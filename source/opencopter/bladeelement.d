module opencopter.bladeelement;

import opencopter.aircraft;
import opencopter.atmosphere;
import opencopter.inflow;
import opencopter.liftmodels;
import opencopter.math;
import opencopter.memory;
import opencopter.trim;
import opencopter.wake;

import std.algorithm;
import std.array;
import std.conv;
import std.math;

extern (C++) void compute_blade_properties(BG, BS, RG, RIS, RS, AS, I, W)(auto ref BG blade, auto ref BS blade_state, auto ref RG rotor, auto ref RIS rotor_input, auto ref RS rotor_state, auto ref AS ac_state, I inflow, auto ref W wake, double time, double dt, size_t rotor_idx, size_t blade_idx)
	if(is_blade_geometry!BG && is_blade_state!BS && is_rotor_geometry!RG && is_rotor_input_state!RIS && is_rotor_state!RS && is_aircraft_state!AS && is_wake!W)
{
	version(LDC) pragma(inline, true);
	version(GNU) pragma(inline, true);

	static import std.math;

	import std.stdio : writeln;

	//immutable cos_azimuth = std.math.cos(blade_state.azimuth);
	//immutable sin_azimuth = std.math.sin(blade_state.azimuth);

	immutable cos_beta = std.math.cos(rotor_input.blade_flapping[blade_idx]);
	immutable sin_beta = std.math.sin(rotor_input.blade_flapping[blade_idx]);

	immutable cos_alpha = rotor_input.cos_aoa;
	immutable sin_alpha = rotor_input.sin_aoa;

	foreach(chunk_idx; 0..blade.chunks.length) {

		immutable Chunk effective_azimuth = blade_state.azimuth - std.math.sgn(rotor_input.angular_velocity)*blade.chunks[chunk_idx].sweep[];
		immutable Chunk cos_sweep = cos(blade.chunks[chunk_idx].sweep);
		immutable Chunk cos_azimuth = cos(effective_azimuth);
		immutable Chunk sin_azimuth = sin(effective_azimuth);

		auto wake_velocities = wake.compute_wake_induced_velocities(blade_state.chunks[chunk_idx].x, blade_state.chunks[chunk_idx].y, blade_state.chunks[chunk_idx].z, ac_state, std.math.abs(rotor_input.angular_velocity), rotor_idx, blade_idx);

		immutable Chunk wake_z = (wake_velocities.v_x[])*sin_alpha + wake_velocities.v_z[]*cos_alpha;
		immutable Chunk wake_y = wake_velocities.v_y[];
		immutable Chunk wake_x = wake_velocities.v_x[]*cos_alpha - wake_velocities.v_z[]*sin_alpha;

		Chunk u_p = -wake_z[] + rotor_state.axial_advance_ratio;
		blade_state.chunks[chunk_idx].u_p[] = u_p[];

		immutable Chunk mu_sin_azimuth = -rotor_state.advance_ratio*sin_azimuth[];
		//immutable Chunk wake_u_t = sin_azimuth[]*wake_x[] + cos_azimuth[]*wake_y[];
		immutable Chunk wake_u_t = cos_azimuth[]*wake_x[] + sin_azimuth[]*wake_y[];
		immutable Chunk sweep_corrected_r = blade.chunks[chunk_idx].r[]*cos_sweep[];
		immutable Chunk u_t = sweep_corrected_r[] + std.math.sgn(rotor_input.angular_velocity)*mu_sin_azimuth[];// + std.math.sgn(rotor_input.angular_velocity)*wake_u_t[]*sin_sweep[];

		Chunk corrected_u_t = u_t[];//.map!(a => a < 0 ? 0 : a).staticArray!Chunk;
		immutable Chunk inflow_angle = atan2(u_p, corrected_u_t);

		blade_state.chunks[chunk_idx].u_t[] = u_t[];
		immutable Chunk plunging_correction = ((rotor_input.blade_flapping_rate[blade_idx]/rotor_input.angular_velocity)*blade.chunks[chunk_idx].r[])/u_t[];
		immutable Chunk theta = (rotor_input.blade_pitches[blade_idx] + blade.chunks[chunk_idx].twist[])[]*cos_sweep[];
		blade_state.chunks[chunk_idx].inflow_angle[] = inflow_angle[];
		blade_state.chunks[chunk_idx].aoa[] = theta[] - inflow_angle[] + plunging_correction[];

		immutable Chunk u_squared = (corrected_u_t[]*corrected_u_t[] + u_p[]*u_p[]);
		immutable Chunk u_inf = sqrt(u_squared);

		immutable Chunk rescaled_u_t = u_inf[]*blade.blade_length/rotor.radius;

		blade_state.circulation_model.compute_bound_circulation_band(blade_state, rescaled_u_t, chunk_idx, std.math.sgn(rotor_input.angular_velocity));

		blade_state.chunks[chunk_idx].aoa_eff[] = -std.math.sgn(rotor_input.angular_velocity)*blade_state.chunks[chunk_idx].gamma[];
		blade_state.chunks[chunk_idx].aoa_eff[] /= (u_inf[]*blade.chunks[chunk_idx].chord[]*PI*2.0*PI);

		auto af_coefficients = blade.airfoil.compute_coeffiecients(chunk_idx, blade_state.chunks[chunk_idx].aoa_eff, zero);

		immutable Chunk dC_L = steady_lift_model(u_p, u_t, af_coefficients.C_l, blade.chunks[chunk_idx].chord)[];

		blade_state.chunks[chunk_idx].dC_L_dot = (dC_L[] - blade_state.chunks[chunk_idx].dC_L[])/dt;
		blade_state.chunks[chunk_idx].dC_L[] = dC_L[];

		immutable Chunk cos_inflow = cos(inflow_angle);
		immutable Chunk cos_collective = std.math.cos(rotor_input.blade_pitches[blade_idx]);
		immutable Chunk sin_collective = std.math.sin(rotor_input.blade_pitches[blade_idx]);
		
		immutable Chunk dC_N = blade_state.chunks[chunk_idx].dC_L[]*cos_collective[];
		immutable Chunk dC_c = -blade_state.chunks[chunk_idx].dC_L[]*sin_collective[];

		immutable Chunk dC_T = blade_state.chunks[chunk_idx].dC_L[]*cos_inflow[];
		blade_state.chunks[chunk_idx].dC_T_dot = (dC_T[] - blade_state.chunks[chunk_idx].dC_T[])/dt;
		blade_state.chunks[chunk_idx].dC_T[] = dC_T[];
		blade_state.chunks[chunk_idx].dC_N[] = dC_N[];
		blade_state.chunks[chunk_idx].dC_c[] = dC_c[];

		blade_state.chunks[chunk_idx].dC_Mx[] = dC_T[]*blade.chunks[chunk_idx].r[]*sin_azimuth[];
		blade_state.chunks[chunk_idx].dC_My[] = dC_T[]*blade.chunks[chunk_idx].r[]*cos_azimuth[];
	}

	//blade_state.set_state_array!"dC_T"([-1.07296e-06, -3.87126e-06, -7.48887e-06, -1.115e-05, -1.42369e-05, -1.58938e-05, -1.49358e-05, -9.98399e-06, 4.18314e-07, 1.77348e-05, 4.33815e-05, 7.86773e-05, 0.000124797, 0.000182727, 0.00025323, 0.00033681, 0.000433695, 0.000543816, 0.000666808, 0.000802011, 0.000948482, 0.00110502, 0.00127019, 0.00144238, 0.00161979, 0.00180053, 0.00198266, 0.0021642, 0.0023432, 0.0025178, 0.00268623, 0.00284688, 0.00299829, 0.00313919, 0.0032685, 0.00338532, 0.00348893, 0.00357878, 0.00365451, 0.00371609, 0.00376475, 0.00380634, 0.00385909, 0.00396032, 0.00414499, 0.0043899, 0.00461532, 0.00476076]);
	blade_state.C_T = integrate_trapaziodal!"dC_T"(blade_state, blade);
	blade_state.C_Mx = integrate_trapaziodal!"dC_Mx"(blade_state, blade);
	blade_state.C_My = integrate_trapaziodal!"dC_My"(blade_state, blade);
}

/++
 +	With a given rotor angual velocity and angular acceleration, compute the lift, torque, power of the rotor.
 +	This is intended to by wrapped in some sort of trim algo.
 +/
extern (C++) void compute_rotor_properties(RG, RS, RIS, AS, I, W)(auto ref RG rotor, auto ref RS rotor_state, auto ref RIS rotor_input, auto ref AS ac_state, I inflow, auto ref W wake, double C_Ti, double time, double dt, size_t rotor_idx)
	if(is_rotor_geometry!RG && is_rotor_input_state!RIS && is_rotor_state!RS && is_aircraft_state!AS && is_wake!W)
{
	version(LDC) pragma(inline, true);
	version(GNU) pragma(inline, true);
	
	rotor_state.C_T = 0.0;
	rotor_state.C_Mx = 0.0;
	rotor_state.C_My = 0.0;

	import std.math : cos, sin, abs;
	import std.stdio : writeln;

	rotor_state.advance_ratio = rotor_input.freestream_velocity*rotor_input.cos_aoa/abs(rotor_input.angular_velocity*rotor.radius);
	rotor_state.axial_advance_ratio = rotor_input.freestream_velocity*rotor_input.sin_aoa/abs(rotor_input.angular_velocity*rotor.radius);

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
			inflow,
			wake,
			time,
			dt,
			rotor_idx,
			blade_idx
		);

		rotor_state.C_T += rotor_state.blade_states[blade_idx].C_T;
		rotor_state.C_Mx += rotor_state.blade_states[blade_idx].C_Mx;
		rotor_state.C_My += rotor_state.blade_states[blade_idx].C_My;
	}
}

void step(I, ArrayContainer AC = ArrayContainer.None)(ref AircraftStateT!AC ac_state, AircraftT!AC aircraft, ref AircraftInputStateT!AC ac_input_state, I[] inflows, ref WakeHistoryT!AC wake_history, immutable Atmosphere atmo, size_t iteration, double dt) {
	
	import std.conv : to;
	import std.math : PI, cos, sin;
	import std.numeric : findRoot;
	import std.stdio : writeln;

	immutable time = iteration.to!double*dt;
	
	foreach(rotor_idx; 0..aircraft.rotors.length) {
		ac_input_state.rotor_inputs[rotor_idx].cos_aoa = cos(ac_input_state.rotor_inputs[rotor_idx].angle_of_attack);
		ac_input_state.rotor_inputs[rotor_idx].sin_aoa = sin(ac_input_state.rotor_inputs[rotor_idx].angle_of_attack);
		
		static import std.math;

		immutable omega_sgn = std.math.sgn(ac_input_state.rotor_inputs[rotor_idx].angular_velocity);

		foreach(blade_idx, ref blade; aircraft.rotors[rotor_idx].blades) {
			ac_state.rotor_states[rotor_idx].blade_states[blade_idx].azimuth = ac_input_state.rotor_inputs[rotor_idx].azimuth + aircraft.rotors[rotor_idx].blades[blade_idx].azimuth_offset;

			immutable cos_beta = std.math.cos(ac_input_state.rotor_inputs[rotor_idx].blade_flapping[blade_idx]);
			immutable sin_beta = std.math.sin(ac_input_state.rotor_inputs[rotor_idx].blade_flapping[blade_idx]);

			immutable cos_azimuth = std.math.cos(ac_state.rotor_states[rotor_idx].blade_states[blade_idx].azimuth);
			immutable sin_azimuth = std.math.sin(ac_state.rotor_states[rotor_idx].blade_states[blade_idx].azimuth);

			immutable cos_alpha = ac_input_state.rotor_inputs[rotor_idx].cos_aoa;
			immutable sin_alpha = ac_input_state.rotor_inputs[rotor_idx].sin_aoa;

			foreach(chunk_idx, ref state_chunk; ac_state.rotor_states[rotor_idx].blade_states[blade_idx].chunks) {
				immutable Chunk x_f = cos_beta*blade.chunks[chunk_idx].r[];
				immutable Chunk z_f = -sin_beta*blade.chunks[chunk_idx].r[];

				immutable Chunk x_tpp = aircraft.rotors[rotor_idx].origin[0] + x_f[]*cos_azimuth + omega_sgn*blade.chunks[chunk_idx].xi[]*sin_azimuth;
				state_chunk.y[] = aircraft.rotors[rotor_idx].origin[1] + x_f[]*sin_azimuth - omega_sgn*blade.chunks[chunk_idx].xi[]*std.math.cos(ac_input_state.rotor_inputs[rotor_idx].blade_pitches[blade_idx])*cos_azimuth;
				immutable Chunk z_tpp = aircraft.rotors[rotor_idx].origin[2] + z_f[] - blade.chunks[chunk_idx].xi[]*std.math.sin(ac_input_state.rotor_inputs[rotor_idx].blade_pitches[blade_idx]);

				state_chunk.x[] = x_tpp[]*cos_alpha + z_tpp[]*sin_alpha;
				state_chunk.z[] = -x_tpp[]*sin_alpha + z_tpp[]*cos_alpha;
			}
		}
	}

	foreach(rotor_idx; 0..aircraft.rotors.length) {
		aircraft.rotors[rotor_idx].compute_rotor_properties(
			ac_state.rotor_states[rotor_idx],
			ac_input_state.rotor_inputs[rotor_idx],
			ac_state,
			inflows[rotor_idx],
			wake_history[0],
			ac_state.rotor_states[rotor_idx].C_T,
			time,
			dt,
			rotor_idx
		);
	}

	aircraft.update_wake(ac_state, ac_input_state, wake_history, inflows, atmo, iteration, dt);

	foreach(rotor_idx; 0..aircraft.rotors.length) {
		inflows[rotor_idx].update(ac_state.rotor_states[rotor_idx].C_T, ac_input_state.rotor_inputs[rotor_idx], ac_state.rotor_states[rotor_idx], ac_state.rotor_states[rotor_idx].advance_ratio, ac_state.rotor_states[rotor_idx].axial_advance_ratio, dt);
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
