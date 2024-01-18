module opencopter.bladeelement;

import opencopter.aircraft;
import opencopter.atmosphere;
import opencopter.inflow;
import opencopter.liftmodels;
import opencopter.math;
import opencopter.memory;
import opencopter.trim;
import opencopter.wake;

import numd.linearalgebra.matrix;

import std.algorithm;
import std.array;
import std.conv;
import std.math;

void compute_blade_properties(BG, BS, RG, RIS, RS, AS, I, WIS, WG, W)(auto ref BG blade, auto ref BS blade_state, auto ref RG rotor, auto ref RIS rotor_input, auto ref RS rotor_state, auto ref AS ac_state, I[] wing_inflows, auto ref WIS wing_input_states, auto ref WG wings, auto ref W wake, double time, double dt, size_t rotor_idx, size_t blade_idx, immutable Atmosphere atmo)
	if(is_blade_geometry!BG && is_blade_state!BS && is_rotor_geometry!RG && is_rotor_input_state!RIS && is_rotor_state!RS && is_aircraft_state!AS && is_wake!W)
{
	version(LDC) pragma(inline, true);
	version(GNU) pragma(inline, true);

	static import std.math;

	import std.stdio : writeln;

	immutable cos_azimuth = std.math.cos(blade_state.azimuth);
	immutable sin_azimuth = std.math.sin(blade_state.azimuth);

	immutable cos_beta = std.math.cos(rotor_input.blade_flapping[blade_idx]);
	immutable sin_beta = std.math.sin(rotor_input.blade_flapping[blade_idx]);

	immutable cos_alpha = rotor_input.cos_aoa;
	immutable sin_alpha = rotor_input.sin_aoa;

	immutable Chunk chunk_of_zeros = 0.0;

	Chunk plunging_correction;

	Chunk wing_v_z = 0.0;
	Chunk wing_v_y = 0.0;
	Chunk wing_v_x = 0.0;

	foreach(chunk_idx; 0..blade.chunks.length) {

		immutable Chunk effective_azimuth = blade_state.azimuth - std.math.sgn(rotor_input.angular_velocity)*blade.chunks[chunk_idx].sweep[];
		immutable Chunk cos_sweep = cos(blade.chunks[chunk_idx].sweep);
		immutable Chunk cos_azimuth = cos(effective_azimuth);
		immutable Chunk sin_azimuth = sin(effective_azimuth);

		auto wake_velocities = wake.compute_wake_induced_velocities(blade_state.chunks[chunk_idx].x, blade_state.chunks[chunk_idx].y, blade_state.chunks[chunk_idx].z, ac_state, rotor_idx, blade_idx, false, false, true);

		auto shed_wake_velocities = wake.compute_wake_induced_velocities(blade_state.chunks[chunk_idx].x, blade_state.chunks[chunk_idx].y, blade_state.chunks[chunk_idx].z, ac_state, rotor_idx, blade_idx, false, true, false);

		auto wake_global_vel = Vector!(4, Chunk)(0);
		auto shed_wake_global_vel = Vector!(4, Chunk)(0);

		wake_global_vel[0][] = wake_velocities.v_x[];
		wake_global_vel[1][] = wake_velocities.v_y[];
		wake_global_vel[2][] = wake_velocities.v_z[];

		shed_wake_global_vel[0][] = shed_wake_velocities.v_x[];
		shed_wake_global_vel[1][] = shed_wake_velocities.v_y[];
		shed_wake_global_vel[2][] = shed_wake_velocities.v_z[];

		immutable Chunk sweep_corrected_r = blade.chunks[chunk_idx].r[]*cos_sweep[];

		immutable total_vel_vec = blade.frame.global_matrix.transpose() * (wake_global_vel + ac_state.freestream);
		immutable shed_vel_vec = blade.frame.global_matrix.transpose() * shed_wake_global_vel;

		blade_state.chunks[chunk_idx].blade_local_vel = total_vel_vec;

		immutable vel_dot = total_vel_vec.dot(blade.chunks[chunk_idx].af_norm);
		immutable projected_vel = total_vel_vec - blade.chunks[chunk_idx].af_norm*vel_dot;

		immutable shed_vel_dot = shed_vel_vec.dot(blade.chunks[chunk_idx].af_norm);
		immutable shed_projected_vel = shed_vel_vec - blade.chunks[chunk_idx].af_norm*shed_vel_dot;

		blade_state.chunks[chunk_idx].projected_vel = projected_vel;


		foreach(wing_idx, wing_inflow; wing_inflows){

			immutable Chunk x_inflow = -x[]*wing_input_states[wing_idx].cos_aoa - z[]*wing_input_states[wing_idx].sin_aoa + wings[wing_idx].origin[0];
			immutable Chunk y_inflow = -y[] + wings[wing_idx].origin[1];
			immutable Chunk z_inflow = -x[]*wing_input_states[wing_idx].sin_aoa + z[]*wing_input_states[wing_idx].cos_aoa + wings[wing_idx].origin[2];

			auto wing_ind_vel = wing_inflow.compute_wing_induced_vel_on_blade(x_inflow,y_inflow,z_inflow);
			
			//immutable aoa = wing_input_states[wing_idx].angle_of_attack;

			immutable w_cos_aoa = wing_input_states[wing_idx].cos_aoa;
			immutable w_sin_aoa = wing_input_states[wing_idx].sin_aoa; 

			wing_v_z[] += (wing_ind_vel.v_z[]*w_cos_aoa - wing_ind_vel.v_x[]*w_sin_aoa)/std.math.abs(rotor_input.angular_velocity); 
			wing_v_x[] += (-wing_ind_vel.v_z[]*w_sin_aoa - wing_ind_vel.v_x[]*w_cos_aoa)/std.math.abs(rotor_input.angular_velocity);
			wing_v_y[] += -wing_ind_vel.v_y[]/std.math.abs(rotor_input.angular_velocity); 
		}

		wake_z[] += wing_v_z[];
		wake_y[] += wing_v_y[];
		wake_x[] += wing_v_x[];

		immutable Chunk wake_z = -projected_vel[2][] + wing_v_z[];

		immutable Chunk u_p = (wake_z[] - shed_projected_vel[2][])/(rotor.radius*abs(rotor_input.angular_velocity));
		blade_state.chunks[chunk_idx].shed_u_p[] = -shed_projected_vel[2][]/(rotor.radius*abs(rotor_input.angular_velocity));
		blade_state.chunks[chunk_idx].u_p[] = u_p[];

		immutable Chunk mu_sin_azimuth = -rotor_state.advance_ratio*sin_azimuth[];
		immutable Chunk u_t = sweep_corrected_r[] + std.math.sgn(rotor_input.angular_velocity)*mu_sin_azimuth[];
		
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

		auto gamma = blade_state.circulation_model.compute_bound_circulation_band(blade_state, chunk_idx, rotor_input.angular_velocity, blade.airfoil.lift_curve_slope(chunk_idx), blade.airfoil.zero_lift_aoa(chunk_idx));

		// Denormalize gamma
		gamma[] *= 0.5 * blade.blade_length * dimensional_u_inf[];

		blade_state.chunks[chunk_idx].d_gamma[] = blade_state.chunks[chunk_idx].gamma[] - gamma[];
		blade_state.chunks[chunk_idx].gamma[] = gamma[];

		blade_state.chunks[chunk_idx].aoa_eff[] = -2.0*std.math.sgn(rotor_input.angular_velocity)*blade_state.chunks[chunk_idx].gamma[];
		blade_state.chunks[chunk_idx].aoa_eff[] /= (u_inf[]*blade.airfoil.lift_curve_slope(chunk_idx)[]*blade.chunks[chunk_idx].chord[]*std.math.abs(rotor_input.angular_velocity)*rotor.radius*rotor.radius);
		blade_state.chunks[chunk_idx].aoa_eff[] += blade.airfoil.zero_lift_aoa(chunk_idx)[];

		auto af_coefficients = blade.airfoil.compute_coeffiecients(chunk_idx, blade_state.chunks[chunk_idx].aoa_eff, M_inf);

		gamma[] = -std.math.sgn(rotor_input.angular_velocity)*0.5*dimensional_u_inf[]*blade.chunks[chunk_idx].chord[]*rotor.radius*af_coefficients.C_l[];

		immutable Chunk dC_L = steady_sectional_model(u_p, u_t, af_coefficients.C_l, blade.chunks[chunk_idx].chord)[];
		immutable Chunk dC_D = steady_sectional_model(u_p, u_t, af_coefficients.C_d, blade.chunks[chunk_idx].chord)[];

		blade_state.chunks[chunk_idx].dC_l[] = af_coefficients.C_l[]*cos(blade_state.chunks[chunk_idx].aoa)[];
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

		immutable Chunk dC_T = (blade_state.chunks[chunk_idx].dC_L[]*cos_inflow[]-blade_state.chunks[chunk_idx].dC_D[]*sin_inflow[]);
		immutable Chunk dC_Db = blade_state.chunks[chunk_idx].dC_L[]*sin_inflow[]+blade_state.chunks[chunk_idx].dC_D[]*cos_inflow[];

		blade_state.chunks[chunk_idx].dC_T_dot = (dC_T[] - blade_state.chunks[chunk_idx].dC_T[])/dt;
		blade_state.chunks[chunk_idx].dC_T[] = dC_T[];
		blade_state.chunks[chunk_idx].dC_Db[] = dC_Db[];
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
}

/++
 +	With a given rotor angual velocity and angular acceleration, compute the lift, torque, power of the rotor.
 +	This is intended to by wrapped in some sort of trim algo.
 +/
void compute_rotor_properties(RG, RS, RIS, AS, WIS, WG, W, I)(auto ref RG rotor, auto ref RS rotor_state, auto ref RIS rotor_input, auto ref AS ac_state,  auto ref WIS wing_input_states, auto ref WG wings, auto ref W wake, auto ref I[] wing_inflows, double C_Ti,double C_Qi, double time, double dt, size_t rotor_idx, immutable Atmosphere atmo)
	if(is_rotor_geometry!RG && is_rotor_input_state!RIS && is_rotor_state!RS && is_aircraft_state!AS && is_wake!W)
{
	version(LDC) pragma(inline, true);
	version(GNU) pragma(inline, true);
	
	rotor_state.C_T = 0.0;
	rotor_state.C_Q = 0.0;
	rotor_state.C_Mx = 0.0;
	rotor_state.C_My = 0.0;

	import std.math : cos, sin, abs;
	import std.stdio : writeln;

	foreach(blade_idx; 0..rotor.blades.length) {
		rotor_state.blade_states[blade_idx].azimuth = rotor_input.azimuth + rotor.blades[blade_idx].azimuth_offset;
	}

	foreach(blade_idx; 0..rotor.blades.length) {
		rotor.blades[blade_idx].compute_blade_properties(
			rotor_state.blade_states[blade_idx],
			rotor,
			rotor_input,
			rotor_state,
			ac_state,
			wing_inflows,
			wing_input_states,
			wings,
			wake,
			time,
			dt,
			rotor_idx,
			blade_idx,
			atmo
		);

		auto blade_frame_forces = Vec4(0, 0, rotor_state.blade_states[blade_idx].C_T, 0);
		auto blade_frame_moments = Vec4(0.0, rotor_state.blade_states[blade_idx].C_My, rotor_state.blade_states[blade_idx].C_Mz, 0);
		
		auto global_frame_forces = rotor.blades[blade_idx].frame.global_matrix*blade_frame_forces;
		auto rotor_frame_forces = rotor.frame.parent.global_matrix.inverse.get()*global_frame_forces;

		auto global_frame_moments = rotor.blades[blade_idx].frame.global_matrix*blade_frame_moments;
		auto rotor_frame_moments = rotor.frame.parent.global_matrix.inverse.get()*global_frame_moments;

		rotor_state.C_T += rotor_frame_forces[2];
		rotor_state.C_Q += rotor_frame_moments[2];
		rotor_state.C_Mx += rotor_frame_moments[0];
		rotor_state.C_My += rotor_frame_moments[1];
	}
}

void print_frame(F)(F frame, int depth = 0) {
	import std.stdio : writeln;
	import std.range : repeat;

	writeln("\t".repeat(depth).join, " ", frame.name, ": ", frame.global_matrix);

	foreach(ref child; frame.children) {
		print_frame(child, depth + 1);
	}
}

void step(ArrayContainer AC = ArrayContainer.None, I)(ref AircraftStateT!AC ac_state, ref AircraftT!AC aircraft, ref AircraftInputStateT!AC ac_input_state, I[] inflows, ref WakeHistoryT!AC wake_history, immutable Atmosphere atmo, size_t iteration, double dt) {
	
	import opencopter.config : chunk_size;

	import std.conv : to;
	import std.math : PI, cos, sin;
	import std.numeric : findRoot;
	import std.stdio : writeln;

	immutable time = iteration.to!double*dt;
	
	aircraft.root_frame.update(Mat4.identity);

	size_t num_rotors = aircraft.rotors.length;

	auto wing_inflows = inflows[num_rotors..$];

	foreach(rotor_idx; 0..aircraft.rotors.length) {

		auto rotor_local_freestream = aircraft.rotors[rotor_idx].frame.parent.global_matrix.inverse.get() * ac_state.freestream;
		
		ac_state.rotor_states[rotor_idx].advance_ratio = abs(rotor_local_freestream[0])/abs(ac_input_state.rotor_inputs[rotor_idx].angular_velocity*aircraft.rotors[rotor_idx].radius);
		ac_state.rotor_states[rotor_idx].axial_advance_ratio = rotor_local_freestream[2]/abs(ac_input_state.rotor_inputs[rotor_idx].angular_velocity*aircraft.rotors[rotor_idx].radius);
		
		foreach(blade_idx, ref blade; aircraft.rotors[rotor_idx].blades) {

			auto local_blade_pos = Vector!(4, Chunk)(0);
			
			foreach(chunk_idx, ref state_chunk; ac_state.rotor_states[rotor_idx].blade_states[blade_idx].chunks) {

				immutable Chunk adjusted_r = blade.chunks[chunk_idx].r[] - blade.r_c;
				local_blade_pos[0][] = adjusted_r[]*aircraft.rotors[rotor_idx].radius;
				local_blade_pos[1][] = blade.chunks[chunk_idx].xi[]*aircraft.rotors[rotor_idx].radius;
				local_blade_pos[1][] = 0;
				local_blade_pos[3][] = 1;

				auto global_blade_pos = blade.frame.global_matrix * local_blade_pos;

				state_chunk.x[] = global_blade_pos[0][];
				state_chunk.y[] = global_blade_pos[1][];
				state_chunk.z[] = global_blade_pos[2][];
			}
		}
	}

	foreach(wing_idx; 0..aircraft.wings.length){
		ac_input_state.wing_inputs[wing_idx].cos_aoa = cos(ac_input_state.wing_inputs[wing_idx].angle_of_attack);
		ac_input_state.wing_inputs[wing_idx].sin_aoa = sin(ac_input_state.wing_inputs[wing_idx].angle_of_attack);
	}

	foreach(wing_idx; 0..aircraft.wings.length){
		ac_input_state.wing_inputs[wing_idx].cos_aoa = cos(ac_input_state.wing_inputs[wing_idx].angle_of_attack);
		ac_input_state.wing_inputs[wing_idx].sin_aoa = sin(ac_input_state.wing_inputs[wing_idx].angle_of_attack);
	}

	foreach(rotor_idx; 0..aircraft.rotors.length) {
		aircraft.rotors[rotor_idx].compute_rotor_properties(
			ac_state.rotor_states[rotor_idx],
			ac_input_state.rotor_inputs[rotor_idx],
			ac_state,			
			ac_input_state.wing_inputs,
			aircraft.wings,
			wake_history[0],
			wing_inflows,
			ac_state.rotor_states[rotor_idx].C_T,
			ac_state.rotor_states[rotor_idx].C_Q,
			time,
			dt,
			rotor_idx,
			atmo
		);
	}

	aircraft.update_wake(ac_state, ac_input_state, wake_history, inflows, atmo, iteration, dt);

	get_ind_vel_on_rotor(ac_state.rotor_states, aircraft.rotors, ac_input_state.rotor_inputs, aircraft.wings, ac_input_state.wing_inputs, ac_state.wing_states, inflows);

	foreach(inflow_idx; num_rotors..inflows.length){
			inflows[inflow_idx].update(ac_input_state, aircraft, inflows,ac_input_state.wing_inputs[0].freestream_velocity, 0.0, 0.0, dt);	
	}

	foreach(rotor_idx; 0..inflows.length) {
		if(rotor_idx >= aircraft.rotors.length) {
			inflows[rotor_idx].update_wing_circulation();
			inflows[rotor_idx].update_wing_dC_L();
		} else {
			//writeln("adv_ratio= ", ac_state.rotor_states[rotor_idx].advance_ratio, "axial_adv_ratio= ",ac_state.rotor_states[rotor_idx].axial_advance_ratio);
			inflows[rotor_idx].update(ac_input_state, aircraft, inflows, ac_input_state.rotor_inputs[0].freestream_velocity, ac_state.rotor_states[rotor_idx].advance_ratio, ac_state.rotor_states[rotor_idx].axial_advance_ratio, dt);
		}
	}
}

import core.thread;
import core.sync.barrier;
//import core.sync.barrier;

auto run(I, ArrayContainer AC = ArrayContainer.None)(ref AircraftStateT!AC ac_state, AircraftT!AC aircraft, ref AircraftInputStateT!AC ac_input_state, I[] inflows, ref WakeHistoryT!AC wake_history, immutable Atmosphere atmo, size_t iteration, double dt, size_t iterations) {
	immutable size_t cores_per_blade = 8;
}

auto process_blade(ArrayContainer AC = ArrayContainer.None)() {

}
