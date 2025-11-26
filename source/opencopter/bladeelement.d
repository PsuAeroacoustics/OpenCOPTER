module opencopter.bladeelement;

import opencopter.aircraft;
import opencopter.atmosphere;
import opencopter.inflow;
import opencopter.liftmodels;
import opencopter.math;
import opencopter.memory;
import opencopter.trim;
import opencopter.wake : WakeHistoryT, WakeT, is_wake, update_wake, compute_wake_induced_velocities, InducedVelocities;
import opencopter.bwi;

import numd.linearalgebra.matrix;

import std.algorithm;
import std.array;
import std.conv;
import std.math;
import std.traits;
import std.typecons;

import core.memory;

struct DynamicInflowWake(RS) {

	RS[] rotor_states;
}

InducedVelocities compute_wake_induced_velocities(W, AS)(auto ref W wake, immutable Chunk x, immutable Chunk y, immutable Chunk z, auto ref AS ac_state, size_t rotor_idx, size_t blade_idx, size_t blade_chunk_idx = 0, bool single_rotor = false, bool shed_only = false, bool tip_only = false, bool trackBWIevents = false)
	if(isInstanceOf!(DynamicInflowWake, W))
{
	import std.stdio : writeln;
	InducedVelocities ret;
	
	ret.v_x[] = 0.0;
	ret.v_y[] = 0.0;
	ret.v_z[] = 0.0;
	
	if(!shed_only) {
		auto global_infow = Vector!(4, Chunk)(0);
		foreach(i_rotor_idx, ref interacting_rotor; wake.rotor_states) {
			auto xyz_chunk = Vector!(4, Chunk)(1);

			xyz_chunk.mData[0][] = x[];
			xyz_chunk.mData[1][] = y[];
			xyz_chunk.mData[2][] = z[];

			auto xyz_tpp = interacting_rotor.inflow_model.frame.inverse_global_matrix * xyz_chunk;

			immutable Chunk lambda_i = interacting_rotor.inflow_model.inflow_at(xyz_tpp)[];

			auto local_inflow = Vector!(4, Chunk)(0);

			local_inflow[2][] = lambda_i[];
			local_inflow[3][] = 0.0;
			global_infow += interacting_rotor.inflow_model.frame.global_matrix * local_inflow;
		}

		ret.v_x[] = global_infow[0][];
		ret.v_y[] = global_infow[1][];
		ret.v_z[] = global_infow[2][];
	}

	return ret;
}
//extern (C++) void compute_blade_properties(BG, BS, RG, RIS, RS, AS, I, W)(auto ref BG blade, auto ref BS blade_state, auto ref RG rotor, auto ref RIS rotor_input, auto ref RS rotor_state, auto ref AS ac_state, I inflow, auto ref W wake, double time, double dt, size_t rotor_idx, size_t blade_idx, immutable Atmosphere atmo, bool trackBWIevents, bool converged)
void compute_blade_properties(BG, BS, RG, RIS, RS, AS, W)   (auto ref BG blade, auto ref BS blade_state, auto ref RG rotor, auto ref RIS rotor_input, auto ref RS rotor_state, auto ref AS ac_state, auto ref W wake, double dt, size_t rotor_idx, size_t blade_idx, immutable Atmosphere atmo, size_t iteration, bool trackBWIevents, bool converged)
	if(is_blade_geometry!BG && is_blade_state!BS && is_rotor_geometry!RG && is_rotor_input_state!RIS && is_rotor_state!RS && is_aircraft_state!AS && is_wake!W)
{
	version(LDC) pragma(inline, true);
	version(GNU) pragma(inline, true);

	static import std.math;

	import std.stdio : writeln;

	debug writeln("SOS: ", atmo.speed_of_sound);
	debug writeln("density: ", atmo.density);
	debug writeln("dynamic_viscosity: ", atmo.dynamic_viscosity);

	debug writeln("rotor state frame name: ", rotor_state.inflow_model.frame.name);

	foreach(chunk_idx; 0..blade.chunks.length) {

		immutable Chunk effective_azimuth = blade_state.azimuth - std.math.sgn(rotor_input.angular_velocity)*blade.chunks[chunk_idx].sweep[];
		immutable Chunk cos_sweep = cos(blade.chunks[chunk_idx].sweep);
		immutable Chunk cos_azimuth = cos(effective_azimuth);
		immutable Chunk sin_azimuth = sin(effective_azimuth);

		auto wake_velocities = wake.compute_wake_induced_velocities(blade_state.chunks[chunk_idx].x, blade_state.chunks[chunk_idx].y, blade_state.chunks[chunk_idx].z, ac_state, rotor_idx, blade_idx, chunk_idx, false, false, true, trackBWIevents);

		auto shed_wake_velocities = wake.compute_wake_induced_velocities(blade_state.chunks[chunk_idx].x, blade_state.chunks[chunk_idx].y, blade_state.chunks[chunk_idx].z, ac_state, rotor_idx, blade_idx, 0, false, true, false, false);

		auto wing_global_infow = Vector!(4, Chunk)(0);

		foreach(ref wing_state ; ac_state.wing_states) {

			auto xyz_chunk = Vector!(4, Chunk)(1);
			xyz_chunk.mData[0][] = blade_state.chunks[chunk_idx].x[];
			xyz_chunk.mData[1][] = blade_state.chunks[chunk_idx].y[];
			xyz_chunk.mData[2][] = blade_state.chunks[chunk_idx].z[];
			debug writeln("/n assigned xyz_chunk: x = ",blade_state.chunks[chunk_idx].x[]);

			debug writeln("wing state inflow frame name: ", wing_state.inflow_model.frame);
			
			auto xyz_tpp = wing_state.inflow_model.frame.inverse_global_matrix * xyz_chunk;
			debug writeln("converted xyz to wing coordinates:", xyz_tpp);

			auto wing_ind_vel = wing_state.inflow_model.compute_wing_induced_vel_on_blade(xyz_tpp[0], xyz_tpp[1], xyz_tpp[2]);
			debug writeln("wing induced velocity calculated, V = ", wing_ind_vel);

			auto local_wing_inflow = Vector!(4, Chunk)(0);
			local_wing_inflow[0][] = wing_ind_vel.v_x[];
			local_wing_inflow[1][] = wing_ind_vel.v_y[];
			local_wing_inflow[2][] = wing_ind_vel.v_z[];

			wing_global_infow += wing_state.inflow_model.frame.global_matrix * local_wing_inflow;
			debug writeln("wing inflow in the global coordinates is calculated to be,", wing_global_infow);
		}

		auto wake_global_vel = Vector!(4, Chunk)(0);
		auto shed_wake_global_vel = Vector!(4, Chunk)(0);

		wake_global_vel[0][] = wake_velocities.v_x[] + wing_global_infow[0][];
		wake_global_vel[1][] = wake_velocities.v_y[] + wing_global_infow[1][];
		wake_global_vel[2][] = wake_velocities.v_z[] + wing_global_infow[2][];

		/*writeln("wake_velocities: ");
		writeln("\tv_x = ", wake_velocities.v_x[]);
		writeln("\tv_y = ", wake_velocities.v_y[]);
		writeln("\tv_z = ", wake_velocities.v_z[]);
		writeln("wing global inflow: ");
		writeln("\tv_x = ", wing_global_infow[0][]);
		writeln("\tv_y = ", wing_global_infow[1][]);
		writeln("\tv_z = ", wing_global_infow[2][]);*/

		shed_wake_global_vel[0][] = shed_wake_velocities.v_x[];
		shed_wake_global_vel[1][] = shed_wake_velocities.v_y[];
		shed_wake_global_vel[2][] = shed_wake_velocities.v_z[];

		/*writeln("shed wake velocities: ");
		writeln("\tv_x = ", shed_wake_velocities.v_x[]);
		writeln("\tv_y = ", shed_wake_velocities.v_x[]);
		writeln("\tv_z = ", shed_wake_velocities.v_x[]);*/

		//immutable total_vel_vec = blade.frame.global_matrix.inverse.get() * (wake_global_vel + ac_state.freestream);
		//immutable shed_vel_vec = blade.frame.global_matrix.inverse.get() * shed_wake_global_vel;
		immutable total_vel_vec = blade.frame.inverse_global_matrix * (wake_global_vel + ac_state.freestream);
		immutable shed_vel_vec = blade.frame.inverse_global_matrix * shed_wake_global_vel;

		blade_state.chunks[chunk_idx].blade_local_vel = total_vel_vec;

		immutable vel_dot = total_vel_vec.dot(blade.chunks[chunk_idx].af_norm);
		immutable projected_vel = total_vel_vec - blade.chunks[chunk_idx].af_norm*vel_dot;

		immutable shed_vel_dot = shed_vel_vec.dot(blade.chunks[chunk_idx].af_norm);
		immutable shed_projected_vel = shed_vel_vec - blade.chunks[chunk_idx].af_norm*shed_vel_dot;

		blade_state.chunks[chunk_idx].projected_vel = /+blade.frame.global_matrix*+/projected_vel;

		immutable Chunk wake_z = -projected_vel[2][];

		immutable Chunk u_p = (wake_z[] - shed_projected_vel[2][])/(rotor.radius*abs(rotor_input.angular_velocity));//*cos_sweep[];

		blade_state.chunks[chunk_idx].shed_u_p[] = -shed_projected_vel[2][]/(rotor.radius*abs(rotor_input.angular_velocity));
		blade_state.chunks[chunk_idx].u_p[] = u_p[];

		immutable Chunk mu_sin_azimuth = -rotor_state.advance_ratio*sin_azimuth[];

		immutable Chunk u_t = (blade.chunks[chunk_idx].r[] + std.math.sgn(rotor_input.angular_velocity)*mu_sin_azimuth[])*cos_sweep[];
		
		immutable Chunk inflow_angle = atan2(u_p, u_t);

		blade_state.chunks[chunk_idx].u_t[] = u_t[];
		immutable Chunk plunging_correction = ((rotor_input.blade_flapping_rate[blade_idx]/abs(rotor_input.angular_velocity))*blade.chunks[chunk_idx].r[])/u_t[];

		immutable Chunk theta = (rotor_input.blade_pitches[blade_idx] + blade.chunks[chunk_idx].twist[])[]*cos_sweep[];
		blade_state.chunks[chunk_idx].theta[] = theta[];
		blade_state.chunks[chunk_idx].inflow_angle[] = inflow_angle[];
		blade_state.chunks[chunk_idx].aoa[] = theta[] - inflow_angle[] - plunging_correction[];

		immutable Chunk u_squared = (u_t[]*u_t[] + u_p[]*u_p[]);
		immutable Chunk u_inf = sqrt(u_squared);
		immutable Chunk dimensional_u_inf = u_inf[] * rotor.radius * abs(rotor_input.angular_velocity);
		immutable Chunk M_inf = dimensional_u_inf[]/atmo.speed_of_sound;

		//writeln("blade sec u_squared = ", u_squared);
		// auto gamma = blade_state.circulation_model.compute_bound_circulation_band(blade_state, chunk_idx, rotor_input.angular_velocity, blade.airfoil.lift_curve_slope(chunk_idx), blade.airfoil.zero_lift_aoa(chunk_idx));

		// //writeln("balde gamma = ", gamma);
		// // Denormalize gamma
		// gamma[] *= 0.5 * blade.blade_length * dimensional_u_inf[];
		// Nitya: Blade circulation dimensionalized here!!

		// blade_state.chunks[chunk_idx].aoa_eff[] = -2.0*std.math.sgn(rotor_input.angular_velocity)*gamma[];

		// blade_state.chunks[chunk_idx].aoa_eff[] /= (u_inf[]*blade.airfoil.lift_curve_slope(chunk_idx)[]*blade.chunks[chunk_idx].chord[]*std.math.abs(rotor_input.angular_velocity)*rotor.radius*rotor.radius);
		// blade_state.chunks[chunk_idx].aoa_eff[] += blade.airfoil.zero_lift_aoa(chunk_idx)[];

		// auto af_coefficients = blade.airfoil.compute_coeffiecients(chunk_idx, blade_state.chunks[chunk_idx].aoa, M_inf);

		// // blade_state.chunks[chunk_idx].d_gamma[] = blade_state.chunks[chunk_idx].gamma[] - gamma[];
		// // blade_state.chunks[chunk_idx].gamma[] = gamma[];

		// immutable Chunk dC_L = steady_sectional_model(u_p, u_t, af_coefficients.C_l, blade.chunks[chunk_idx].chord)[];
		// immutable Chunk dC_D = steady_sectional_model(u_p, u_t, af_coefficients.C_d, blade.chunks[chunk_idx].chord)[];

		// blade_state.chunks[chunk_idx].dC_l[] = af_coefficients.C_l[];
		// blade_state.chunks[chunk_idx].dC_d[] = af_coefficients.C_d[];

		// blade_state.chunks[chunk_idx].dC_L_dot = (dC_L[] - blade_state.chunks[chunk_idx].dC_L[])/dt;
		// blade_state.chunks[chunk_idx].dC_L[] = dC_L[];
		// blade_state.chunks[chunk_idx].dC_D[] = dC_D[];

		// immutable Chunk cos_inflow = cos(inflow_angle);
		// immutable Chunk sin_inflow = sin(inflow_angle);

		// immutable Chunk cos_collective = std.math.cos(rotor_input.blade_pitches[blade_idx]);
		// immutable Chunk sin_collective = std.math.sin(rotor_input.blade_pitches[blade_idx]);

		// immutable Chunk dC_N = blade_state.chunks[chunk_idx].dC_L[]*cos_collective[];
		// immutable Chunk dC_c = -blade_state.chunks[chunk_idx].dC_L[]*sin_collective[];

		// immutable Chunk dC_T = (blade_state.chunks[chunk_idx].dC_L[]*cos_inflow[] - blade_state.chunks[chunk_idx].dC_D[]*sin_inflow[]);
		// immutable Chunk dC_Db = blade_state.chunks[chunk_idx].dC_L[]*sin_inflow[] + blade_state.chunks[chunk_idx].dC_D[]*cos_inflow[];

		// blade_state.chunks[chunk_idx].dC_T_dot = (dC_T[] - blade_state.chunks[chunk_idx].dC_T[])/dt;
		// blade_state.chunks[chunk_idx].dC_T[] = dC_T[];
		// blade_state.chunks[chunk_idx].dT[] = dC_T[] * atmo.density*PI*rotor.radius^^4.0*abs(rotor_input.angular_velocity)^^2.0;
		// blade_state.chunks[chunk_idx].dC_Db[] = dC_Db[];
		// blade_state.chunks[chunk_idx].dC_Db_induced[] = blade_state.chunks[chunk_idx].dC_L[]*sin_inflow[];
		// blade_state.chunks[chunk_idx].dC_Db_profile[] = blade_state.chunks[chunk_idx].dC_D[]*cos_inflow[];
		// blade_state.chunks[chunk_idx].dC_N[] = dC_N[];
		// blade_state.chunks[chunk_idx].dC_c[] = dC_c[];

		// blade_state.chunks[chunk_idx].dC_My[] = dC_T[]*blade.chunks[chunk_idx].r[];
		// blade_state.chunks[chunk_idx].dC_Mz[] = dC_Db[]*blade.chunks[chunk_idx].r[];
	}

	foreach(chunk_idx; 0..blade.chunks.length) {
		auto gamma = blade_state.circulation_model.compute_bound_circulation_band(blade_state, chunk_idx, rotor_input.angular_velocity, blade.airfoil.lift_curve_slope(chunk_idx), blade.airfoil.zero_lift_aoa(chunk_idx));

		immutable Chunk u_p = blade_state.chunks[chunk_idx].u_p[];
		immutable Chunk u_t = blade_state.chunks[chunk_idx].u_t[];

		immutable Chunk u_squared = (u_t[]*u_t[] + u_p[]*u_p[]);
		immutable Chunk u_inf = sqrt(u_squared);
		immutable Chunk dimensional_u_inf = u_inf[] * rotor.radius * abs(rotor_input.angular_velocity);
		immutable Chunk M_inf = dimensional_u_inf[]/atmo.speed_of_sound;
		//writeln("balde gamma = ", gamma);
		// Denormalize gamma
		gamma[] *= 0.5 * blade.blade_length * dimensional_u_inf[];// * abs(rotor_input.angular_velocity);

		blade_state.chunks[chunk_idx].d_gamma[] = blade_state.chunks[chunk_idx].gamma[] - gamma[];
		blade_state.chunks[chunk_idx].gamma[] = gamma[];

		blade_state.chunks[chunk_idx].aoa_eff[] = 2.0*std.math.sgn(rotor_input.angular_velocity)*gamma[];

		blade_state.chunks[chunk_idx].aoa_eff[] /= (u_inf[]*blade.airfoil.lift_curve_slope(chunk_idx)[]*blade.chunks[chunk_idx].chord[]*std.math.abs(rotor_input.angular_velocity)*rotor.radius*rotor.radius);
		blade_state.chunks[chunk_idx].aoa_eff[] += blade.airfoil.zero_lift_aoa(chunk_idx)[];

		auto af_coefficients = blade.airfoil.compute_coeffiecients(chunk_idx, blade_state.chunks[chunk_idx].aoa_eff, M_inf);

		// blade_state.chunks[chunk_idx].d_gamma[] = blade_state.chunks[chunk_idx].gamma[] - gamma[];
		// blade_state.chunks[chunk_idx].gamma[] = gamma[];

		immutable Chunk dC_L = steady_sectional_model(u_p, u_t, af_coefficients.C_l, blade.chunks[chunk_idx].chord)[];
		immutable Chunk dC_D = steady_sectional_model(u_p, u_t, af_coefficients.C_d, blade.chunks[chunk_idx].chord)[];

		immutable Chunk plunging_correction = ((rotor_input.blade_flapping_rate[blade_idx]/abs(rotor_input.angular_velocity))*blade.chunks[chunk_idx].r[])/u_t[];
		immutable Chunk inflow_angle = blade_state.chunks[chunk_idx].theta[] - blade_state.chunks[chunk_idx].aoa_eff[] - plunging_correction[];

		blade_state.chunks[chunk_idx].effective_inflow_angle[] = inflow_angle[];

		blade_state.chunks[chunk_idx].dC_l[] = af_coefficients.C_l[];
		blade_state.chunks[chunk_idx].dC_d[] = af_coefficients.C_d[];

		blade_state.chunks[chunk_idx].dC_L_dot = (dC_L[] - blade_state.chunks[chunk_idx].dC_L[])/dt;
		blade_state.chunks[chunk_idx].dC_L[] = dC_L[];
		blade_state.chunks[chunk_idx].dC_D[] = dC_D[];

		immutable Chunk cos_inflow = cos(inflow_angle);
		immutable Chunk sin_inflow = sin(inflow_angle);

		immutable Chunk cos_collective = std.math.cos(rotor_input.blade_pitches[blade_idx]);
		immutable Chunk sin_collective = std.math.sin(rotor_input.blade_pitches[blade_idx]);

		immutable Chunk dC_N = blade_state.chunks[chunk_idx].dC_L[]*cos_collective[];
		immutable Chunk dC_c = -blade_state.chunks[chunk_idx].dC_L[]*sin_collective[];

		immutable Chunk dC_T = (blade_state.chunks[chunk_idx].dC_L[]*cos_inflow[] - blade_state.chunks[chunk_idx].dC_D[]*sin_inflow[]);
		immutable Chunk dC_Db = blade_state.chunks[chunk_idx].dC_L[]*sin_inflow[] + blade_state.chunks[chunk_idx].dC_D[]*cos_inflow[];

		blade_state.chunks[chunk_idx].dC_T_dot = (dC_T[] - blade_state.chunks[chunk_idx].dC_T[])/dt;
		blade_state.chunks[chunk_idx].dC_T[] = dC_T[];
		blade_state.chunks[chunk_idx].dT[] = dC_T[] * atmo.density*PI*rotor.radius^^4.0*abs(rotor_input.angular_velocity)^^2.0;
		blade_state.chunks[chunk_idx].dC_Db[] = dC_Db[];
		blade_state.chunks[chunk_idx].dC_Db_induced[] = blade_state.chunks[chunk_idx].dC_L[]*sin_inflow[];
		blade_state.chunks[chunk_idx].dC_Db_profile[] = blade_state.chunks[chunk_idx].dC_D[]*cos_inflow[];
		blade_state.chunks[chunk_idx].dC_N[] = dC_N[];
		blade_state.chunks[chunk_idx].dC_c[] = dC_c[];

		blade_state.chunks[chunk_idx].dC_My[] = dC_T[]*blade.chunks[chunk_idx].r[];
		blade_state.chunks[chunk_idx].dC_Mz[] = dC_Db[]*blade.chunks[chunk_idx].r[];
	}

	blade_state.chunks[$-1].d_gamma[$-1] = 0;
	blade_state.chunks[0].d_gamma[0] = 0;

	blade_state.C_T = integrate_trapaziodal!"dC_T"(blade_state, blade);
	blade_state.C_Mz = integrate_trapaziodal!"dC_Mz"(blade_state, blade);
	blade_state.C_My = integrate_trapaziodal!"dC_My"(blade_state, blade);

	// Nitya, 09.14
	
	static if(!isInstanceOf!(DynamicInflowWake, W))
	{
		if(converged && trackBWIevents){
			auto normal = Vector!(4, Chunk)(0);
			normal[2][] = 1;
			auto updated_normal = blade.frame.global_matrix * normal;
			double[3] normalVec;
			normalVec[0] = updated_normal[0][0];
			normalVec[1] = updated_normal[1][0];
			normalVec[2] = updated_normal[2][0];
			//double[] r = get_geometry_array!"r"(blade);
			//writeln("1. bladeElement, before  iteration:", iteration);
			// foreach (i_rotor_idx; 0..wake.rotor_wakes.length){
			// 	//writeln("wake.rotor_wakes[i_rotor_idx].interaction_perRotor[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction.length:", wake.rotor_wakes[i_rotor_idx].interaction_perRotor[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction.length);
			// 	foreach (i_blade_idx; 0..wake.rotor_wakes[i_rotor_idx].interaction_perRotor[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction.length) {
			// 		/*if(wake.rotor_wakes[i_rotor_idx].interaction_perRotor[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction[i_blade_idx].interaction_pts.empty){
			// 			writeln("no interaction point");
			// 		} else{
			// 			wake.rotor_wakes[i_rotor_idx].interaction_perRotor[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction[i_blade_idx].interaction_pts.clear();
			// 		}*/
			// 		wake.rotor_wakes[i_rotor_idx].interaction_perRotor[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction[i_blade_idx].interaction_pts.clear();
			// 	} 
			// }
			calculate_BWI_points(wake, blade_state, rotor_idx, blade_idx, blade, normalVec, iteration);
			// writeln("1. bladeElement, after iteration:", iteration);
		}
	}
}

/++
 +	With a given rotor angual velocity and angular acceleration, compute the lift, torque, power of the rotor.
 +	This is intended to by wrapped in some sort of trim algo.
 +/
 //extern (C++) void compute_rotor_properties(RG, RS, RIS, AS, I, W)      (auto ref RG rotor, auto ref RS rotor_state, auto ref RIS rotor_input, auto ref AS ac_state, I inflow, auto ref W wake, double C_Ti,double C_Qi, double time, double dt, size_t rotor_idx, immutable Atmosphere atmo, bool trackBWIevents, bool converged)
void compute_rotor_properties(RG, RS, RIS, AS, WIS, WG, W)(auto ref RG rotor, auto ref RS rotor_state, auto ref RIS rotor_input, auto ref AS ac_state,  auto ref WIS wing_input_states, auto ref WG wings, auto ref W wake, double C_Ti,double C_Qi, size_t iteration, double dt, size_t rotor_idx, immutable Atmosphere atmo, bool trackBWIevents, bool converged)
	if(is_rotor_geometry!RG && is_rotor_input_state!RIS && is_rotor_state!RS && is_aircraft_state!AS && is_wake!W)
{
	version(LDC) pragma(inline, true);
	version(GNU) pragma(inline, true);

	double C_T = 0;
	rotor_state.C_Q = 0.0;
	rotor_state.C_Mx = 0.0;
	rotor_state.C_My = 0.0;

	import std.math : cos, sin, abs;
	import std.stdio : writeln;

	//debug writeln("wing geometry frame: ", wings[0].frame);
	//debug writeln("wing geometry frame: ", wings[0].frame.name);

	foreach(blade_idx; 0..rotor.blades.length) {
		rotor_state.blade_states[blade_idx].azimuth = rotor_input.azimuth + rotor.blades[blade_idx].azimuth_offset;
		// Nitya: blade azimuth calculated here! 
	}

	auto dynWake = DynamicInflowWake!(PointerTarget!RS)(ac_state.rotor_states.data);

	Chunk[] backup_CT = new Chunk[rotor.blades[0].chunks.length];
	Chunk[] backup_T = new Chunk[rotor.blades[0].chunks.length];
	// Chunk[] backup_CMy = new Chunk[rotor.blades[0].chunks.length];
	// Chunk[] backup_CMz = new Chunk[rotor.blades[0].chunks.length];
	Chunk[] backup_up = new Chunk[rotor.blades[0].chunks.length];
	Chunk[] backup_inflow_angle = new Chunk[rotor.blades[0].chunks.length];
	Chunk[] backup_gamma = new Chunk[rotor.blades[0].chunks.length];
	Chunk[] backup_d_gamma = new Chunk[rotor.blades[0].chunks.length];

	foreach(blade_idx; 0..rotor.blades.length) {

		// if(iteration > 0) {
		// 	//debug writeln("1. iteration:", iteration);
		// 	//writeln("\n blade properites for blade ", blade_idx);
		// 	rotor.blades[blade_idx].compute_blade_properties(
		// 		rotor_state.blade_states[blade_idx],
		// 		rotor,
		// 		rotor_input,
		// 		rotor_state,
		// 		ac_state,
		// 		dynWake,
		// 		dt,
		// 		rotor_idx,
		// 		blade_idx,
		// 		atmo,
		// 		iteration,
		// 		trackBWIevents,
		// 		converged
		// 	);
		// 	//writeln("blade state calculation done \n");

		// 	auto blade_frame_forces = Vec4(0, 0, rotor_state.blade_states[blade_idx].C_T, 0);
		// 	auto blade_frame_moments = Vec4(0.0, rotor_state.blade_states[blade_idx].C_My, rotor_state.blade_states[blade_idx].C_Mz, 0.0);
			
		// 	auto global_frame_forces = rotor.blades[blade_idx].frame.global_matrix*blade_frame_forces;
		// 	auto rotor_frame_forces = rotor.frame.parent.global_matrix.inverse.get()*global_frame_forces;

		// 	auto global_frame_moments = rotor.blades[blade_idx].frame.global_matrix*blade_frame_moments;
		// 	auto rotor_frame_moments = rotor.frame.parent.global_matrix.inverse.get()*global_frame_moments;

		// 	C_T += rotor_frame_forces[2];
		// 	rotor_state.C_Q += rotor_frame_moments[2];

			
		// 	rotor_state.C_Mx += rotor_frame_moments[0];
		// 	rotor_state.C_My += rotor_frame_moments[1];
		// }

		// foreach(chunk_idx, ref blade_chunk; rotor_state.blade_states[blade_idx].chunks) {
		// 	backup_CT[chunk_idx][] = blade_chunk.dC_T[];
		// 	backup_T[chunk_idx][] = blade_chunk.dT[];
		// 	// backup_CMy[chunk_idx][] = blade_chunk.dC_My[];
		// 	// backup_CMz[chunk_idx][] = blade_chunk.dC_Mz[];
		// 	backup_up[chunk_idx][] = blade_chunk.u_p[];
		// 	backup_inflow_angle[chunk_idx][] = blade_chunk.inflow_angle[];
		// 	backup_gamma[chunk_idx][] = blade_chunk.gamma[];
		// 	backup_d_gamma[chunk_idx][] = blade_chunk.d_gamma[];
		// 	blade_chunk.dynamic_u_p[] = blade_chunk.u_p[];
		// 	blade_chunk.dynamic_dC_Db_induced[] = blade_chunk.dC_Db_induced[];
		// 	blade_chunk.dynamic_dC_Db_profile[] = blade_chunk.dC_Db_profile[];
		// 	blade_chunk.dynamic_u_p[] = blade_chunk.u_p[];
		// 	blade_chunk.dynamic_inflow_angle[] = blade_chunk.inflow_angle[];
		// 	blade_chunk.aoa_eff[] = blade_chunk.aoa[];
		// }
		//debug writeln("2. iteration:", iteration);
		//debug writeln("blade_idx:", blade_idx, "rotor_idx:", rotor_idx);
		//writeln("\n blade properites for blade ", blade_idx);
		rotor.blades[blade_idx].compute_blade_properties(
			rotor_state.blade_states[blade_idx],
			rotor,
			rotor_input,
			rotor_state,
			ac_state,
			wake,
			dt,
			rotor_idx,
			blade_idx,
			atmo, 
			iteration,
			trackBWIevents,
			converged
		);
		//writeln("blade state calculation done \n");

		// foreach(chunk_idx, ref blade_chunk; rotor_state.blade_states[blade_idx].chunks) {
		// 	if(iteration > 0) {
		// 		blade_chunk.dC_T[] = backup_CT[chunk_idx][];
		// 		blade_chunk.dT[] = backup_T[chunk_idx][];
		// 		blade_chunk.gamma[] = backup_gamma[chunk_idx][];
		// 		//blade_chunk.u_p[] = backup_up[chunk_idx][];
		// 		//blade_chunk.inflow_angle[] = backup_inflow_angle[chunk_idx][];
		// 		blade_chunk.d_gamma[] = backup_d_gamma[chunk_idx][];
		// 	}
		// 	// blade_chunk.dC_My[] = backup_CMy[chunk_idx][];
		// 	// blade_chunk.dC_Mz[] = backup_CMz[chunk_idx][];
		// }

		// if(iteration == 0) {
			auto blade_frame_forces = Vec4(0, 0, rotor_state.blade_states[blade_idx].C_T, 0);
			auto blade_frame_moments = Vec4(0.0, rotor_state.blade_states[blade_idx].C_My, rotor_state.blade_states[blade_idx].C_Mz, 0.0);
			
			auto global_frame_forces = rotor.blades[blade_idx].frame.global_matrix*blade_frame_forces;
			//auto rotor_frame_forces = rotor.frame.parent.global_matrix.inverse.get()*global_frame_forces;
			auto rotor_frame_forces = rotor.frame.parent.inverse_global_matrix*global_frame_forces;

			auto global_frame_moments = rotor.blades[blade_idx].frame.global_matrix*blade_frame_moments;
			//auto rotor_frame_moments = rotor.frame.parent.global_matrix.inverse.get()*global_frame_moments;
			auto rotor_frame_moments = rotor.frame.parent.inverse_global_matrix*global_frame_moments;

			C_T += rotor_frame_forces[2];
			rotor_state.C_Q += rotor_frame_moments[2];

			rotor_state.C_Mx += rotor_frame_moments[0];
			rotor_state.C_My += rotor_frame_moments[1];
		// }
	}

	rotor_state.C_T = C_T;

	//writeln("rotor_C_T = ", rotor_state.C_T);
}

void print_frame(F)(F frame, int depth = 0) {
	import std.stdio : writeln;
	import std.range : repeat;

	//writeln("\t".repeat(depth).join, " ", frame.name, ": ", frame.local_matrix);
	writeln("\t".repeat(depth).join, " ", frame.name, ": ", frame.global_matrix);
	writeln("\t".repeat(depth).join, " ", frame.name, ": ", frame.inverse_global_matrix);

	foreach(ref child; frame.children) {
		print_frame(child, depth + 1);
	}
}

//void step(I, ArrayContainer AC = ArrayContainer.None)(ref AircraftStateT!AC ac_state, AircraftT!AC aircraft, ref AircraftInputStateT!AC ac_input_state, I[] inflows, ref WakeHistoryT!AC wake_history, immutable Atmosphere atmo, size_t iteration, double dt, bool trackBWIevents, bool converged) {
void step(ArrayContainer AC = ArrayContainer.None)(ref AircraftStateT!AC ac_state, ref AircraftT!AC aircraft, ref AircraftInputStateT!AC ac_input_state, ref WakeHistoryT!AC wake_history, immutable Atmosphere atmo, size_t iteration, double dt, bool trackBWIevents, bool converged) {
	
	import opencopter.config : chunk_size;

	import std.conv : to;
	import std.math : PI, cos, sin;
	import std.numeric : findRoot;
	import std.stdio : writeln;

	//GC.collect();
	
	aircraft.root_frame.update(Mat4.identity);
	//aircraft.root_frame.print_frame;

	foreach(rotor_idx; 0..aircraft.rotors.length) {

		//auto rotor_local_freestream = aircraft.rotors[rotor_idx].frame.parent.global_matrix.inverse.get() * ac_state.freestream;
		auto rotor_local_freestream = aircraft.rotors[rotor_idx].frame.parent.inverse_global_matrix * ac_state.freestream;
		
		ac_state.rotor_states[rotor_idx].advance_ratio = abs(rotor_local_freestream[0])/abs(ac_input_state.rotor_inputs[rotor_idx].angular_velocity*aircraft.rotors[rotor_idx].radius);
		ac_state.rotor_states[rotor_idx].axial_advance_ratio = rotor_local_freestream[2]/abs(ac_input_state.rotor_inputs[rotor_idx].angular_velocity*aircraft.rotors[rotor_idx].radius);
		
		foreach(blade_idx, ref blade; aircraft.rotors[rotor_idx].blades) {

			auto local_blade_pos = Vector!(4, Chunk)(0);
			
			foreach(chunk_idx, ref state_chunk; ac_state.rotor_states[rotor_idx].blade_states[blade_idx].chunks) {

				immutable Chunk adjusted_r = blade.chunks[chunk_idx].r[] - blade.r_c;
				local_blade_pos[0][] = adjusted_r[]*aircraft.rotors[rotor_idx].radius;
				local_blade_pos[1][] = blade.chunks[chunk_idx].xi[]*aircraft.rotors[rotor_idx].radius;
				local_blade_pos[2][] = 0;
				local_blade_pos[3][] = 1;

				auto global_blade_pos = blade.frame.global_matrix * local_blade_pos;

				state_chunk.x[] = global_blade_pos[0][];
				state_chunk.y[] = global_blade_pos[1][];
				state_chunk.z[] = global_blade_pos[2][];
			}
		}
	}

	foreach(rotor_idx; 0..aircraft.rotors.length) {
		//writeln("rotor properties for rotor ", rotor_idx);
		aircraft.rotors[rotor_idx].compute_rotor_properties(
			ac_state.rotor_states[rotor_idx],
			ac_input_state.rotor_inputs[rotor_idx],
			ac_state,
			ac_input_state.wing_inputs,
			aircraft.wings,
			wake_history[0],
			ac_state.rotor_states[rotor_idx].C_T,
			ac_state.rotor_states[rotor_idx].C_Q,
			iteration,
			dt,
			rotor_idx,
			atmo,
			trackBWIevents,
			converged
		);
	}
	
	// Nitya: HOW?? - you can call a function by referencing first argument in the function too!
	//aircraft.update_wake(ac_state, ac_input_state, wake_history, inflows, atmo, iteration, dt);
	aircraft.update_wake(ac_state, ac_input_state, wake_history, atmo, iteration, dt);
	
	debug writeln("aircraft wake updated");
	auto ac_forces = Vec4(0.0);
	//foreach(rotor_idx; 0..aircraft.rotors.length) {
	foreach(rotor_idx, ref rotor_state ; ac_state.rotor_states) {
		//inflows[rotor_idx].update(ac_state.rotor_states[rotor_idx].C_T, ac_input_state.rotor_inputs[rotor_idx], ac_state.rotor_states[rotor_idx], ac_state.rotor_states[rotor_idx].advance_ratio, ac_state.rotor_states[rotor_idx].axial_advance_ratio, &ac_state, dt);
		rotor_state.inflow_model.update(ac_state, dt);
		foreach(blade_idx, ref blade; aircraft.rotors[rotor_idx].blades) {
			foreach(chunk_idx, ref state_chunk; ac_state.rotor_states[rotor_idx].blade_states[blade_idx].chunks) {
				state_chunk.x_old[] = state_chunk.x[];
				state_chunk.y_old[] = state_chunk.y[];
				state_chunk.z_old[] = state_chunk.z[];
			}
		}

		auto rotor_forces = Vec4(0.0);

		rotor_forces[2] = rotor_state.C_T*atmo.density*PI*aircraft.rotors[rotor_idx].radius^^4.0*abs(ac_input_state.rotor_inputs[rotor_idx].angular_velocity)^^2.0;

		auto global_rotor_forces = aircraft.rotors[rotor_idx].frame.global_matrix*rotor_forces;
		//ac_forces += aircraft.root_frame.global_matrix.inverse.get * global_rotor_forces;
		ac_forces += aircraft.root_frame.inverse_global_matrix * global_rotor_forces;
	}

	// debug writeln("\n rotor properties computed");
	// aircraft.update_wake(ac_state, ac_input_state, wake_history, atmo, iteration, dt);


	
	// foreach(r_idx, ref rotor_state ; ac_state.rotor_states) {
	// 	rotor_state.inflow_model.update(ac_state, dt);

	// 	auto rotor_forces = Vec4(0.0);

	// 	rotor_forces[2] = rotor_state.C_T*atmo.density*PI*aircraft.rotors[r_idx].radius^^4.0*abs(ac_input_state.rotor_inputs[r_idx].angular_velocity)^^2.0;

	// 	auto global_rotor_forces = aircraft.rotors[r_idx].frame.global_matrix*rotor_forces;
	// 	ac_forces += aircraft.root_frame.global_matrix.inverse.get * global_rotor_forces;
	// }
	debug writeln("\n rotor inflows updated");
	ac_state.forces = ac_forces;

	foreach(ref wing_state ; ac_state.wing_states) {
		wing_state.inflow_model.update(ac_state, dt);
	}

	debug writeln("\n wing inflows updated");
}
