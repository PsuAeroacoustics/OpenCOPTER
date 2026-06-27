/*
 * OpenCOPTER C++ Example
 *
 * Mirrors examples/python/example.py - single rotor forward flight simulation
 * with Huang-Peters dynamic inflow and free-wake tracking.
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
  const int iterations = 5400;
  const size_t wake_history_length = 1 * 1024;
  const size_t requested_elements = 45;

  const size_t num_rotors = 1;
  const size_t num_blades = 4;
  const double R = 2.0;

  const double theta_75 = 2.0 * (M_PI / 180.0);

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
  std::vector<double> r = oc::generate_radius_points(requested_elements, r_c);
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
  // Build blade geometry and frame hierarchy
  // ============================================================
  oc::Vec3 rotor_origin{0, 0, 0};

  /* Load airfoil polar data.
     AirfoilModel is move-only so we store it in a vector that outlives the
     blades. BladeAirfoil also must outlive all BladeGeometry objects. */
  std::string polar_path = "../polars/NACA23012mod_1000000_polar.dat";
  std::vector<oc::AirfoilModel> airfoil_models;
  airfoil_models.emplace_back(
      oc::AirfoilModel::aero_das_from_xfoil_polar(polar_path, 0.12));

  std::vector<size_t> airfoil_extents{0, 47};
  oc::BladeAirfoil blade_airfoil =
      oc::BladeAirfoil::create(airfoil_models, airfoil_extents);

  /* Average chord for solidity (computed from chord distribution) */
  double avg_chord = R * std::accumulate(c.begin(), c.end(), 0.0) /
                     static_cast<double>(c.size());

  /* Storage for blades and frame tree per blade:
     azimuth_offset_frame -> root_cutout_frame -> blade_frame     */
  std::vector<oc::BladeGeometry> blades;
  blades.reserve(num_blades);

  std::vector<oc::Frame> azimuth_frames(num_blades);
  std::vector<oc::Frame> cutout_frames(num_blades);
  std::vector<oc::Frame> blade_frames(num_blades);

  for (size_t b_idx = 0; b_idx < num_blades; ++b_idx) {
    oc::BladeGeometry blade(static_cast<size_t>(elements), 0.0, avg_chord,
                            blade_airfoil, r_c);
    blade.set_blade_length(R * (1.0 - r_c));

    /* Set geometric properties */
    blade.set_radius(r);
    blade.set_twist(twist);
    blade.set_alpha_0(alpha_0);
    blade.set_C_l_alpha(C_l_alpha);
    blade.set_chord(c);
    blade.set_sweep(sweep);
    blade.set_xi(xi);
    blade.set_xi_p(xi_p);

    blade.compute_vectors();

    /* Build frame hierarchy for this blade */
    std::string az_name = "blade_" + std::to_string(b_idx) + "_azimuth";
    azimuth_frames[b_idx] =
        oc::Frame({0, 0, 1}, static_cast<double>(b_idx) * d_azimuth, {0, 0, 0},
                  az_name, oc::FrameType::Connection);

    std::string cut_name = "blade_" + std::to_string(b_idx) + "_cutout";
    cutout_frames[b_idx] = oc::Frame({1, 0, 0}, 0.0, {R * r_c, 0, 0}, cut_name,
                                     oc::FrameType::Connection);

    std::string bl_name = "blade_" + std::to_string(b_idx);
    blade_frames[b_idx] =
        oc::Frame({1, 0, 0}, 0.0, {0, 0, 0}, bl_name, oc::FrameType::Blade);

    /* Parent-child links: azimuth -> cutout -> blade */
    const oc::Frame *az_children[] = {&cutout_frames[b_idx]};
    azimuth_frames[b_idx].set_children(az_children);

    const oc::Frame *cut_children[] = {&blade_frames[b_idx]};
    cutout_frames[b_idx].set_children(cut_children);

    blade.set_frame(blade_frames[b_idx]);

    blades.push_back(std::move(blade));
  }

  // ============================================================
  // Build rotor geometry
  // ============================================================
  oc::RotorGeometry rotor(num_blades, rotor_origin, R, 0.0);

  double solidity = static_cast<double>(num_blades) * avg_chord / (M_PI * R);
  rotor.set_solidity(solidity);

  /* Assign blades to rotor */
  std::vector<const oc::BladeGeometry *> blade_ptrs;
  for (auto &b : blades) {
    blade_ptrs.push_back(&b);
  }
  rotor.set_blades(blade_ptrs);

  /* Rotor frame - parent of all azimuth offset frames */
  oc::Frame rotor_frame({1, 0, 0}, 0.0, {0, 0, 0}, "rotor_0",
                        oc::FrameType::Rotor);
  std::vector<const oc::Frame *> az_ptrs;
  for (auto &af : azimuth_frames) {
    az_ptrs.push_back(&af);
  }
  rotor_frame.set_children(az_ptrs);
  rotor.set_frame(rotor_frame);

  // ============================================================
  // Aircraft container
  // ============================================================
  oc::Aircraft aircraft(num_rotors, 0);

  oc::Frame root_frame = aircraft.root_frame();
  root_frame.set_name("Aircraft frame");
  root_frame.set_frame_type(oc::FrameType::Aircraft);

  /* Apply angle of attack rotation about Y axis */
  root_frame.rotate({0, 1, 0}, -aoa);
  root_frame.update(oc::mat4_identity());

  /* Assign rotor(s) to aircraft */
  std::vector<const oc::RotorGeometry *> rotor_ptrs{&rotor};
  aircraft.set_rotors(rotor_ptrs);

  /* Link aircraft root -> rotor frame */
  const oc::Frame *root_children[] = {&rotor_frame};
  root_frame.set_children(root_children);

  // ============================================================
  // Input state (kinematics fed into aero model)
  // ============================================================
  oc::AircraftInputState ac_input_state(num_rotors, {num_blades}, 0);

  /* Configure rotor 0 input */
  oc::RotorInputState rotor_input = ac_input_state.get_rotor_input(0);
  rotor_input.set_angular_velocity(omega);
  rotor_input.set_angular_accel(0.0);
  rotor_input.set_azimuth(0.0);

  std::vector<double> r_0_values(num_blades, 1.0 * avg_chord);
  std::vector<double> flapping_values(num_blades, 0.0);
  std::vector<double> flapping_rate_values(num_blades, 0.0);
  rotor_input.set_r_0(r_0_values);
  rotor_input.set_blade_flapping(flapping_values);
  rotor_input.set_blade_flapping_rate(flapping_rate_values);

  /* Set collective pitch per blade */
  for (size_t b_idx = 0; b_idx < num_blades; ++b_idx) {
    ac_input_state.set_blade_pitch(0, b_idx, theta_75);
  }

  std::cout << "r_0 values set for rotor 0" << std::endl;

  // ============================================================
  // Huang-Peters dynamic inflow model
  // ============================================================
  oc::HuangPetersInflow inflow_0(4, 2, rotor, rotor_input, dt);

  /* Build std::span<Inflow*> for the AircraftState constructor.
     The API expects explicit spans, not raw arrays. */
  oc::Inflow *rotor_inflow_arr[] = {&inflow_0};
  std::span<oc::Inflow *> rotor_inflows_span(rotor_inflow_arr);
  /* Empty span for wing inflows (no wings in this example).
     Use default constructor which produces a nullptr/0-size span. */
  std::span<oc::Inflow *> wing_inflows_span{};

  // ============================================================
  // Wake history (free-wake storage)
  // ============================================================
  oc::WakeHistory wake_history(num_rotors, num_blades, wake_history_length, 2,
                               elements, {shed_history}, {1}, 6.5e-5, false);

  // ============================================================
  // Aircraft aerodynamic state
  // ============================================================
  oc::AircraftState ac_state(num_rotors, std::vector<size_t>{num_blades},
                             elements, 0, std::vector<size_t>{}, 0,
                             0, /* no wings */
                             aircraft, rotor_inflows_span, wing_inflows_span,
                             oc::direction_counter_clockwise());

  ac_state.set_freestream(oc::Vec4{V_inf, 0, 0, 0});

  // ============================================================
  // VTK output setup
  // ============================================================
  oc::VtkRotor vtk_rotor_0 = oc::VtkRotor::build(rotor);
  const oc::VtkRotor *vtk_ptrs[] = {&vtk_rotor_0};

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

    /* Advance rotor azimuth by one timestep */
    oc::basic_aircraft_rotor_dynamics(ac_input_state, dt);

    /* Perform simulation step (BLADE element + wake update + inflow advance) */
    oc::simulation_step(ac_state, aircraft, ac_input_state, wake_history, atmo,
                        static_cast<size_t>(iteration), dt, false, false);

    /* Write VTK output for last 360 iterations */
    if (iteration > (iterations - 360)) {
      oc::write_rotors_vtu("rotor", static_cast<size_t>(iteration), vtk_ptrs,
                           ac_state, aircraft);

      oc::Wake wake = wake_history.get_wake(0);
      oc::VtkWake vtk_wake = oc::VtkWake::build(wake);
      oc::write_wake_vtu("wake", static_cast<size_t>(iteration), vtk_wake,
                         wake);
    }
  }

  // ============================================================
  // Results
  // ============================================================
  std::cout << "rotor 0 C_T: " << ac_state.rotor_C_T(0) << std::endl;

  /* Extract and print wake vortex data */
  double max_dim_1 = -std::numeric_limits<double>::infinity();
  double min_dim_1 = std::numeric_limits<double>::infinity();
  double max_dim_2 = -std::numeric_limits<double>::infinity();
  double min_dim_2 = std::numeric_limits<double>::infinity();

  oc::Wake final_wake = wake_history.get_wake(0);
  oc::RotorWake rotor_wake = final_wake.get_rotor_wake(0);

  for (size_t b_idx = 0; b_idx < num_blades; ++b_idx) {
    oc::VortexFilament tip_vortex = rotor_wake.get_tip_vortex(b_idx);

    std::vector<double> x = tip_vortex.x(static_cast<size_t>(elements));
    std::vector<double> z = tip_vortex.z(static_cast<size_t>(elements));

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

  return 0;
}