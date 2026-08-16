/*
 * OpenCOPTER C++ Example
 *
 * Mirrors examples/python/example.py - single rotor forward flight simulation
 * with Huang-Peters dynamic inflow.
 */

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

#include <opencopter>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace oc = opencopter;

int main() {
	// ============================================================
	// Configuration
	// ============================================================
	const int iterations = 361;
	const size_t wake_history_length = 1 * 1024;
	const size_t requested_elements = 45;

	const size_t num_rotors = 1;
	const size_t num_blades = 4;
	const double R = 2.0;

	const double theta_75 = 3.2 * (M_PI / 180.0);

	const double density = 1.125;
	const double omega = 109.12;
	const double AR = 16.5;
	const double sos = 343.0;

	const double V_inf = 32.9;
	const double aoa = 6.5 * (M_PI / 180.0);
	const double theta_tw_1 = -6.5 * (M_PI / 180.0);

	const double d_psi = 1.0;
	const double dt = d_psi * (M_PI / 180.0) / std::abs(omega);

	const size_t shed_history = static_cast<size_t>(std::round(45.0 / d_psi));

	const double d_azimuth = 2.0 * M_PI / num_blades;

	/* Root cutout */
	const double r_c = 0.22;

	// ============================================================
	// Spanwise distributions (normalized by rotor radius)
	// ============================================================
	auto r = oc::generate_radius_points(requested_elements, r_c);
	const size_t elements = r.size();
	std::cout << "requested_elements: " << requested_elements
			<< ", actual elements: " << elements << std::endl;

	const double chord_val = 1.0 / AR;
	std::vector<double> c(elements, chord_val);
	std::vector<double> alpha_0(elements, 0.0);
	std::vector<double> xi(elements, 0.0);
	std::vector<double> xi_p(elements, 0.0);

	std::vector<double> twist(elements);
	for (size_t i = 0; i < elements; ++i) {
		twist[i] = (r[i] - 0.75) * theta_tw_1;
		//std::cout << "r[i]: " << r[i] << " twist[" << i << "]: " << twist[i] << std::endl;
	}

	const double cl_alpha_val = 2.0 * M_PI;
	std::vector<double> C_l_alpha(elements, cl_alpha_val);
	std::vector<double> sweep(elements, 0.0);

	// ============================================================
	// Atmosphere
	// ============================================================
	oc::AtmosphereData atmo{density, 18.03e-6,
							density == 0 ? 0.0 : 18.03e-6 / density, sos};

	std::cout << "Freestream vel: " << V_inf << " m/s" << std::endl;

	// ============================================================
	// Aircraft container - created FIRST so root_frame exists as parent
	// (Mirrors Python: aircraft = Aircraft(num_rotors, 0))
	// ============================================================
	oc::Aircraft aircraft(num_rotors, 0);
	oc::Frame aircraft_root = aircraft.root_frame();
	aircraft_root.set_name("Aircraft frame");
	aircraft_root.set_frame_type(oc::FrameType::Aircraft);

	/* Apply angle of attack rotation about Y axis */
	aircraft_root.rotate({0, 1, 0}, -aoa);
	aircraft_root.update(oc::mat4_identity());

	// ============================================================
	// Load airfoil polar data
	// The airfoil_models vector MUST outlive all BladeGeometry objects.
	// ============================================================
	std::string polar_path = "../polars/NACA23012mod_1000000_polar.dat";
	std::vector<oc::AirfoilModel> airfoil_models;
	airfoil_models.emplace_back(
		oc::AirfoilModel::aero_das_from_xfoil_polar(polar_path, 0.12));

	std::vector<size_t> airfoil_extents{0, 47};
	oc::BladeAirfoil blade_airfoil =
		oc::BladeAirfoil::create(airfoil_models, airfoil_extents);

	/* Average chord for solidity */
	double avg_chord = R * std::accumulate(c.begin(), c.end(), 0.0) /
						static_cast<double>(c.size());

	// ============================================================
	// Build blade geometries
	// ============================================================
	std::vector<oc::BladeGeometry> blades;
	blades.reserve(num_blades);

	for (size_t b_idx = 0; b_idx < num_blades; ++b_idx) {
		oc::BladeGeometry blade(static_cast<size_t>(elements), 0.0, avg_chord,
								blade_airfoil, r_c);
		blade.set_blade_length(R * (1.0 - r_c));

		blade.set_radius(r);
		blade.set_twist(twist);
		blade.set_alpha_0(alpha_0);
		blade.set_C_l_alpha(C_l_alpha);
		blade.set_chord(c);
		blade.set_sweep(sweep);
		blade.set_xi(xi);
		blade.set_xi_p(xi_p);

		blade.compute_vectors();

		blades.push_back(std::move(blade));
	}

	// ============================================================
	// Build frame hierarchy:
	//   aircraft_root -> rotor_frame -> azimuth_frames[b]
	//                                         -> cutout_frames[b]
	//                                                      -> blade_frames[b]
	//
	// NOTE: The C++ Frame API does NOT support passing parent at
	// construction. Parent-child relationships are established via
	// set_children() which sets children AND links the parent pointer.
	// This matches the D code behavior where setting .children also
	// sets each child's .parent = this.
	// ============================================================

	/* IMPORTANT: The D backend pins each Frame in GC.addRoot at creation.
 		When frames are linked into a parent-child hierarchy, destroying an
 		individual child frame (calling GC.removeRoot while the parent still
 		holds a D reference) corrupts the pinned list -> heap crash.

 		Solution: Use owned_=true only for the Aircraft which owns the entire
 		tree. All Frame objects created here are marked non-owning so their
 		destructors do NOT call oc_frame_destroy/GC.removeRoot. The frames
 		live on the stack and are kept alive by the D GC through parent refs. */

	/* Create rotor frame with aircraft_root as parent */
	oc::Frame rotor_fixed_frame({1, 0, 0}, 0.0, {0, 0, 0}, &aircraft_root, "rotor_0_fixed", oc::FrameType::Connection);
	oc::Frame rotor_frame({1, 0, 0}, 0.0, {0, 0, 0}, &rotor_fixed_frame, "rotor_0", oc::FrameType::Rotor);

	auto parent_frame = rotor_frame.parent();
	
	const oc::Frame* rotor_frame_ptr = &rotor_frame;
	rotor_fixed_frame.set_children(std::span<const oc::Frame*>(&rotor_frame_ptr, 1));

	std::vector<oc::Frame> azimuth_frames;
	std::vector<oc::Frame> cutout_frames;
	std::vector<oc::Frame> blade_frames;
	azimuth_frames.reserve(num_blades);
	cutout_frames.reserve(num_blades);
	blade_frames.reserve(num_blades);

	for (size_t b_idx = 0; b_idx < num_blades; ++b_idx) {
		std::string az_name = "blade_" + std::to_string(b_idx) + "_azimuth";
		/* azimuth frame: parent is rotor_frame */
		azimuth_frames.emplace_back(
			oc::Frame({0, 0, 1}, static_cast<double>(b_idx) * d_azimuth, {0, 0, 0},
						&rotor_frame, az_name, oc::FrameType::Connection));

		std::string cut_name = "blade_" + std::to_string(b_idx) + "_cutout";
		/* cutout frame: parent is azimuth_frames[b_idx] */
		cutout_frames.emplace_back(
			oc::Frame({1, 0, 0}, 0.0, {R * r_c, 0, 0}, &azimuth_frames[b_idx], cut_name,
						oc::FrameType::Connection));

		std::string bl_name = "blade_" + std::to_string(b_idx);
		/* blade frame: parent is cutout_frames[b_idx] */
		blade_frames.emplace_back(
			oc::Frame({1, 0, 0}, 0.0, {0, 0, 0}, &cutout_frames[b_idx], bl_name, oc::FrameType::Blade));

		/* Assign the leaf blade frame to the blade geometry */
		blades[b_idx].set_frame(blade_frames[b_idx]);
	}

	/* Link children arrays using set_children() for bidirectional references.
 		In the C API, oc_frame_set_children sets child.parent = this for each child,
 		ensuring both parent->children[] and child->parent are consistent. */

	/* blade_frames[b] are children of cutout_frames[b] */
	for (size_t b_idx = 0; b_idx < num_blades; ++b_idx) {
		const oc::Frame* blade_ptr = &blade_frames[b_idx];
		cutout_frames[b_idx].set_children(std::span<const oc::Frame*>(
			&blade_ptr, 1));
	}

	/* cutout_frames[b] are children of azimuth_frames[b] */
	for (size_t b_idx = 0; b_idx < num_blades; ++b_idx) {
		const oc::Frame* cutout_ptr = &cutout_frames[b_idx];
		azimuth_frames[b_idx].set_children(std::span<const oc::Frame*>(
			&cutout_ptr, 1));
	}

	/* azimuth_frames[] are children of rotor_frame */
	{
		std::vector<const oc::Frame*> az_ptrs;
		for (auto& af : azimuth_frames) az_ptrs.push_back(&af);
		rotor_frame.set_children(az_ptrs);
	}

	// ============================================================
	// Build rotor geometry
	// ============================================================
	oc::RotorGeometry rotor(num_blades, {0, 0, 0}, R, 0.0);

	double solidity = static_cast<double>(num_blades) * avg_chord / (M_PI * R);
	std::cout << "solidity: " << solidity << std::endl;
	rotor.set_solidity(solidity);

	/* Assign blades to rotor */
	{
		std::vector<const oc::BladeGeometry*> blade_ptrs;
		for (auto& b : blades) blade_ptrs.push_back(&b);
		rotor.set_blades(blade_ptrs);
	}

	rotor.set_frame(rotor_frame);

	/* Assign rotor to aircraft - the RotorGeometry must live on the stack
		since the Aircraft stores a pointer (non-owning reference). */
	{
		const oc::RotorGeometry* rp = &rotor;
		aircraft.set_rotors(std::span<const oc::RotorGeometry*>(&rp, 1));
	}

	/* Link aircraft_root -> rotor_frame. After this call the D GC chain is:
		aircraft (root frame) -> rotor_frame -> azimuth_frames -> cutout_frames -> blade_frames */
	{
		const oc::Frame* rf = &rotor_fixed_frame;
		aircraft_root.set_children(std::span<const oc::Frame*>(&rf, 1));
	}

	// ============================================================
	// Input state (kinematics fed into aero model)
	// ============================================================
	oc::AircraftInputState ac_input_state(num_rotors, {num_blades}, 0);

	/* Get rotor input handle and configure it */
	oc::RotorInputState rotor_input = ac_input_state.get_rotor_input(0);
	rotor_input.set_angular_velocity(omega);
	rotor_input.set_angular_accel(0.0);
	rotor_input.set_azimuth(0.0);

	{
		std::vector<double> r_0_vals(num_blades, 1.0 * avg_chord);
		std::cout << "r_0: " << r_0_vals[0] << std::endl;
		std::vector<double> flap_vals(num_blades, 0.0);
		std::vector<double> flap_rate_vals(num_blades, 0.0);
		rotor_input.set_r_0(r_0_vals);
		rotor_input.set_blade_flapping(flap_vals);
		rotor_input.set_blade_flapping_rate(flap_rate_vals);
	}

	for (size_t b_idx = 0; b_idx < num_blades; ++b_idx) {
		ac_input_state.set_blade_pitch(0, b_idx, theta_75);
	}

	std::cout << "r_0 values set for rotor 0" << std::endl;

	// ============================================================
	// Huang-Peters dynamic inflow model
	// ============================================================
	oc::HuangPetersInflow inflow_0(4, 2, rotor, rotor_input, dt);

	/* Build vectors for AircraftState constructor */
	std::vector<oc::Inflow*> rotor_inflows = {&inflow_0};
	std::vector<oc::Inflow*> wing_inflows;

	// ============================================================
	// Wake history (free-wake storage)
	// ============================================================
	oc::WakeHistory wake_history(num_rotors, num_blades, wake_history_length, 2, elements,
						{shed_history}, {1}, 6.5e-5, false);

	// ============================================================
	// Aircraft aerodynamic state
	// ============================================================
	oc::AircraftState ac_state(num_rotors, {num_blades}, elements, 0, {}, 0, 0, aircraft,
						rotor_inflows, wing_inflows,
						oc::direction_counter_clockwise());

	ac_state.set_freestream(oc::Vec4{V_inf, 0, 0, 0});

	// ============================================================
	// VTK output setup (build once before simulation loop)
	// ============================================================
	oc::VtkRotor vtk_rotor_0 = oc::VtkRotor::build(rotor);
	if (!vtk_rotor_0) {
		std::cerr << "WARNING: VtkRotor::build returned null — rotor VTK output will be skipped" << std::endl;
	} else {
		std::cout << "VTK rotor created successfully" << std::endl;
	}

	// Wake VTK — grab the first wake from history to build base structure
	oc::Wake initial_wake = wake_history.get_wake(0);
	oc::VtkWake vtk_wake;
	if (initial_wake) {
		vtk_wake = oc::VtkWake::build(initial_wake);
		if (!vtk_wake) {
			std::cerr << "WARNING: VtkWake::build returned null — wake VTK output will be skipped" << std::endl;
		} else {
			std::cout << "VTK wake created successfully" << std::endl;
		}
	} else {
		std::cerr << "WARNING: wake_history.get_wake(0) returned null — wake VTK output will be skipped" << std::endl;
	}

	// ============================================================
	// Simulation loop
	// ============================================================
	auto start_time = std::chrono::steady_clock::now();

	for (int iteration = 0; iteration < iterations; ++iteration) {
		if (iteration % 100 == 0) {
			auto now = std::chrono::steady_clock::now();
			double elapsed = std::chrono::duration<double>(now - start_time).count();
			std::cout << elapsed << ": iteration: " << iteration << std::endl;
			start_time = now;
		}

		oc::basic_aircraft_rotor_dynamics(ac_input_state, dt);

		rotor_frame.set_rotation(oc::Vec3(0, 0, 1), ac_input_state.get_rotor_input(0).azimuth());

		oc::simulation_step(ac_state, aircraft, ac_input_state, wake_history, atmo,
							static_cast<size_t>(iteration), dt, false, false);

		if (vtk_rotor_0 && iteration > (iterations - 360)) {
			std::vector<oc::VtkRotor> vtk_arr = {vtk_rotor_0};
			oc::write_rotors_vtu("rotor", static_cast<size_t>(iteration), vtk_arr,
								ac_state, aircraft);
		}
		if (vtk_wake && iteration > (iterations - 360)) {
			oc::Wake wake = wake_history.get_wake(0);
			oc::write_wake_vtu("wake", static_cast<size_t>(iteration), vtk_wake, wake);
		}
	}

	// ============================================================
	// Results
	// ============================================================
	std::cout << "rotor 0 C_T: " << ac_state.rotor_C_T(0) << std::endl;

	double max_dim_1 = -std::numeric_limits<double>::infinity();
	double min_dim_1 = std::numeric_limits<double>::infinity();
	double max_dim_2 = -std::numeric_limits<double>::infinity();
	double min_dim_2 = std::numeric_limits<double>::infinity();

	oc::Wake final_wake = wake_history.get_wake(0);
	oc::RotorWake rotor_wake = final_wake.get_rotor_wake(0);

	for (size_t b_idx = 0; b_idx < num_blades; ++b_idx) {
	oc::VortexFilament tip_vortex = rotor_wake.get_tip_vortex(b_idx);

	std::vector<double> x = tip_vortex.x(elements);
	std::vector<double> z = tip_vortex.z(elements);

	double x_max = *std::max_element(x.begin(), x.end());
	double x_min = *std::min_element(x.begin(), x.end());
	double z_max = *std::max_element(z.begin(), z.end());
	double z_min = *std::min_element(z.begin(), z.end());

	max_dim_2 = std::max(max_dim_2, z_max);
	min_dim_2 = std::min(min_dim_2, z_min);
	max_dim_1 = std::max(max_dim_1, x_max);
	min_dim_1 = std::min(min_dim_1, x_min);

	std::cout << "Blade " << b_idx << " tip wake: x=[" << x_min << ", " << x_max
				<< "] z=[" << z_min << ", " << z_max << "]" << std::endl;
	}

	double span_1 = max_dim_1 - min_dim_1;
	double span_2 = max_dim_2 - min_dim_2;
	std::cout << "Total wake bounding box: x=" << span_1 << ", z=" << span_2
			<< std::endl;

	delete r.data();

	return 0;
}