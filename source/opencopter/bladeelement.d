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

void compute_blade_properties(alias lift_model, BG, BS, RG, RIS, RS, AS, WIS, WG, W)(auto ref BG blade, auto ref BS blade_state, auto ref RG rotor, auto ref RIS rotor_input, auto ref RS rotor_state, auto ref AS ac_state, Inflow[] wing_inflows, auto ref WIS[] wing_input_states, auto ref WG[] wings, auto ref W wake, double time, double dt, size_t rotor_idx, size_t blade_idx)
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
		immutable Chunk x_f = cos_beta*blade.chunks[chunk_idx].r[];
		immutable Chunk z_f = -sin_beta*blade.chunks[chunk_idx].r[];

		immutable Chunk x_tpp = rotor.origin[0] + x_f[]*cos_azimuth;
		immutable Chunk y = rotor.origin[1] + x_f[]*sin_azimuth;
		immutable Chunk z_tpp = rotor.origin[2] + z_f[];

		immutable Chunk x = x_tpp[]*cos_alpha + z_tpp[]*sin_alpha;
		immutable Chunk z = -x_tpp[]*sin_alpha + z_tpp[]*cos_alpha;
		
		auto wake_velocities = wake.compute_wake_induced_velocities(x, y, z, ac_state, std.math.abs(rotor_input.angular_velocity), rotor_idx, blade_idx);

		Chunk wake_z = (wake_velocities.v_x[]*sin_alpha + wake_velocities.v_z[]*cos_alpha);
		Chunk wake_y = wake_velocities.v_y[];
		Chunk wake_x = wake_velocities.v_x[]*cos_alpha - wake_velocities.v_z[]*sin_alpha;
		
		//writeln("wake_x = ", wake_x, "wake_y = ", wake_y, "wake_z = ", wake_z);
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

		immutable Chunk u_p = -wake_z[] +  rotor_state.axial_advance_ratio;
		//writeln("wake_z = ", wake_z, "\twing_v_z= ", wing_v_z);
		blade_state.chunks[chunk_idx].u_p[] = u_p[];
		//writeln("wake_z = ", wake_z, "\twing_v_z= ", wing_v_z, "\taxial_adv_ratio= ", rotor_state.axial_advance_ratio);
		immutable mu_sin_azimuth = -rotor_state.advance_ratio*sin_azimuth;
		immutable Chunk wake_u_t = sin_azimuth*wake_x[] + cos_azimuth*wake_y[];
		//writeln("wake_u_t = ", wake_u_t);
		immutable Chunk u_t = blade.chunks[chunk_idx].r[] + std.math.sgn(rotor_input.angular_velocity)*mu_sin_azimuth + /+std.math.sgn(rotor_input.angular_velocity)*+/wake_u_t[];

		immutable Chunk corrected_u_t = u_t[];
		immutable Chunk inflow_angle = atan2(u_p, corrected_u_t);
		//writeln("inflow_angle = ", inflow_angle);
		blade_state.chunks[chunk_idx].u_t[] = u_t[];
		if(u_t[0] == 0.0){
			plunging_correction[] = 0.0; //((rotor_input.blade_flapping_rate[blade_idx]/rotor_input.angular_velocity)*blade.chunks[chunk_idx].r[])/u_t[];
		} else {
			plunging_correction[] = ((rotor_input.blade_flapping_rate[blade_idx]/rotor_input.angular_velocity)*blade.chunks[chunk_idx].r[])/u_t[];	
		}
		//immutable Chunk plunging_correction = ((rotor_input.blade_flapping_rate[blade_idx]/rotor_input.angular_velocity)*blade.chunks[chunk_idx].r[])/u_t[];
		immutable Chunk theta = rotor_input.blade_pitches[blade_idx] + blade.chunks[chunk_idx].twist[];
		blade_state.chunks[chunk_idx].inflow_angle[] = inflow_angle[];
		blade_state.chunks[chunk_idx].aoa[] = theta[] - inflow_angle[] + plunging_correction[];
		//writeln("chunk aoa = ", blade_state.chunks[chunk_idx].aoa);
		//writeln("r = ", blade.chunks[chunk_idx].r[]);
		//writeln("u_t = ",u_t);
		auto af_coefficients = blade.airfoil.compute_coeffiecients(chunk_idx, blade_state.chunks[chunk_idx].aoa, zero);
		//writeln("af_coefficient = ", af_coefficients);
		immutable Chunk rescaled_u_t = u_t[]/blade.average_chord;
		//writeln("u_t = ", u_t);
		blade_state.circulation_model.compute_bound_circulation_band(blade_state, rescaled_u_t, chunk_idx);
		//writeln("u_p= ", u_p, "\tu_t= ", u_t);
		immutable Chunk dC_L = steady_lift_model(u_p, u_t, af_coefficients.C_l, blade.chunks[chunk_idx].chord)[];

		blade_state.chunks[chunk_idx].dC_L_dot = (dC_L[] - blade_state.chunks[chunk_idx].dC_L[])/dt;
		blade_state.chunks[chunk_idx].dC_L[] = dC_L[];

		immutable Chunk cos_inflow = cos(inflow_angle);
		
		immutable Chunk dC_T = blade_state.chunks[chunk_idx].dC_L[]*cos_inflow[];
		blade_state.chunks[chunk_idx].dC_T_dot = (dC_T[] - blade_state.chunks[chunk_idx].dC_T[])/dt;
		blade_state.chunks[chunk_idx].dC_T[] = dC_T[];

		blade_state.chunks[chunk_idx].dC_Mx[] = dC_T[]*blade.chunks[chunk_idx].r[]*sin_azimuth;
		blade_state.chunks[chunk_idx].dC_My[] = dC_T[]*blade.chunks[chunk_idx].r[]*cos_azimuth;
	}

	blade_state.C_T = integrate_trapaziodal!"dC_T"(blade_state, blade);
	blade_state.C_Mx = integrate_trapaziodal!"dC_Mx"(blade_state, blade);
	blade_state.C_My = integrate_trapaziodal!"dC_My"(blade_state, blade);
}

/++
 +	With a given rotor angual velocity and angular acceleration, compute the lift, torque, power of the rotor.
 +	This is intended to by wrapped in some sort of trim algo.
 +/
void compute_rotor_properties(alias lift_model, RG, RS, RIS, AS, WIS, WG, W)(auto ref RG rotor, auto ref RS rotor_state, auto ref RIS rotor_input, auto ref AS ac_state,  auto ref WIS[] wing_input_states, auto ref WG[] wings, auto ref W wake, Inflow[] wing_inflows, double C_Ti, double time, double dt, size_t rotor_idx)
	if(is_rotor_geometry!RG && is_rotor_input_state!RIS && is_rotor_state!RS && is_aircraft_state!AS && is_wake!W)
{
	version(LDC) pragma(inline, true);
	version(GNU) pragma(inline, true);
	
	rotor_state.C_T = 0.0;
	rotor_state.C_Mx = 0.0;
	rotor_state.C_My = 0.0;

	import std.math : cos, sin, abs;
	import std.stdio : writeln;

	//rotor_state.advance_ratio = rotor_input.freestream_velocity*rotor_input.cos_aoa/abs(rotor_input.angular_velocity*rotor.radius);
	//rotor_state.axial_advance_ratio = rotor_input.freestream_velocity*rotor_input.sin_aoa/abs(rotor_input.angular_velocity*rotor.radius);

	/*writeln("freestream velocity:", rotor_input.freestream_velocity);
	writeln("cos_aoa:", rotor_input.cos_aoa);
	writeln("sin_aoa:", rotor_input.sin_aoa);
	writeln("omega:", rotor_input.angular_velocity);
	writeln("radius:", rotor.radius);

	writeln("adv ratio: ", rotor_state.advance_ratio);
	writeln("axial adv ratio: ", rotor_state.axial_advance_ratio);*/


	foreach(blade_idx; 0..rotor.blades.length) {
		rotor_state.blade_states[blade_idx].azimuth = rotor_input.azimuth + rotor.blades[blade_idx].azimuth_offset;
	}

	foreach(blade_idx; 0..rotor.blades.length) {

		rotor.blades[blade_idx].compute_blade_properties!lift_model(
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
			blade_idx
		);

		rotor_state.C_T += rotor_state.blade_states[blade_idx].C_T;
		rotor_state.C_Mx += rotor_state.blade_states[blade_idx].C_Mx;
		rotor_state.C_My += rotor_state.blade_states[blade_idx].C_My;
	}
}

void step(ArrayContainer AC = ArrayContainer.None)(ref AircraftStateT!AC ac_state, AircraftT!AC aircraft, ref AircraftInputStateT!AC ac_input_state, Inflow[] inflows, ref WakeHistoryT!AC wake_history, immutable Atmosphere atmo, size_t iteration, double dt) {
	
	import std.conv : to;
	import std.math : PI, cos, sin;
	import std.numeric : findRoot;
	import std.stdio : writeln;

	immutable time = iteration.to!double*dt;

	size_t num_rotors = aircraft.rotors.length;

	auto wing_inflows = inflows[num_rotors..$];
	
	foreach(rotor_idx; 0..aircraft.rotors.length) {
		ac_input_state.rotor_inputs[rotor_idx].cos_aoa = cos(ac_input_state.rotor_inputs[rotor_idx].angle_of_attack);
		ac_input_state.rotor_inputs[rotor_idx].sin_aoa = sin(ac_input_state.rotor_inputs[rotor_idx].angle_of_attack);
	}

	foreach(wing_idx; 0..aircraft.wings.length){
		ac_input_state.wing_inputs[wing_idx].cos_aoa = cos(ac_input_state.wing_inputs[wing_idx].angle_of_attack);
		ac_input_state.wing_inputs[wing_idx].sin_aoa = sin(ac_input_state.wing_inputs[wing_idx].angle_of_attack);
	}

	foreach(rotor_idx; 0..aircraft.rotors.length) {
		aircraft.rotors[rotor_idx].compute_rotor_properties!steady_lift_model(
			ac_state.rotor_states[rotor_idx],
			ac_input_state.rotor_inputs[rotor_idx],
			ac_state,			
			ac_input_state.wing_inputs,
			aircraft.wings,
			wake_history[0],
			wing_inflows,
			ac_state.rotor_states[rotor_idx].C_T,
			time,
			dt,
			rotor_idx
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
