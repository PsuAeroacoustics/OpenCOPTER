module opencopter.trim;

import opencopter.aircraft;
import opencopter.atmosphere;
import opencopter.bladeelement;
import opencopter.inflow;
import opencopter.liftmodels;
import opencopter.wake;

import std.meta;

/+alias trim_list = aliasSeqOf!(["trim_to_thrust_rpm", "no_trim"]);

alias TrimFunction = extern (C++) @nogc void function(ref AircraftState ac_state, Aircraft aircraft, ref AircraftInputState ac_input_state, ref AircraftInputState last_ac_input_state, ref Atmosphere atmo, Inflow inflow, ref Wake wake, size_t rotor_idx, double time, double dt);

private @nogc double trim_to_thrust_rpm_impl(alias lift_model)(double anglular_velocity, ref AircraftState ac_state, Aircraft aircraft, ref AircraftInputState ac_input_state, ref AircraftInputState last_ac_input_state, ref Atmosphere atmo, Inflow inflow, ref Wake wake, size_t rotor_idx, double time, double dt) {

	import std.conv : to;
	import std.math : PI;
	import std.numeric : findRoot;
	import std.stdio : writeln;


	ac_input_state.rotor_inputs[rotor_idx].angular_velocity = anglular_velocity;

	immutable R = aircraft.rotors[rotor_idx].radius;
	immutable R_squared = R*R;
	immutable Omega = ac_input_state.rotor_inputs[rotor_idx].angular_velocity;

	immutable desired_C_T = ac_input_state.rotor_inputs[rotor_idx].required_thrust/(atmo.density*PI*R_squared*R_squared*Omega*Omega);

	aircraft.rotors[rotor_idx].compute_rotor_properties!lift_model(
		ac_state.rotor_states[rotor_idx],
		ac_input_state.rotor_inputs[rotor_idx],
		ac_state,
		inflow,
		wake,
		desired_C_T,
		time,
		dt
	);

	immutable rotor_thrust = atmo.density*PI*R_squared*R_squared*Omega*Omega*ac_state.rotor_states[rotor_idx].C_T;

	immutable ret = rotor_thrust - ac_input_state.rotor_inputs[rotor_idx].required_thrust;
	//debug writeln("rotor_thrust: ", rotor_thrust);
	//debug writeln("ret: ", ret);
	return rotor_thrust - ac_input_state.rotor_inputs[rotor_idx].required_thrust;
}

extern (C++) @nogc void trim_to_thrust_rpm(alias lift_model)(ref AircraftState ac_state, Aircraft aircraft, ref AircraftInputState ac_input_state, ref AircraftInputState last_ac_input_state, ref Atmosphere atmo, Inflow inflow, ref Wake wake, size_t rotor_idx, double time, double dt) {

	import std.conv : to;
	import std.math : PI;
	import std.numeric : findRoot;
	import std.stdio : writeln;

	@nogc double trim(double anglular_velocity) {
		return trim_to_thrust_rpm_impl!lift_model(anglular_velocity, ac_state, aircraft, ac_input_state, last_ac_input_state, atmo, inflow, wake, rotor_idx, time, dt);
	}

	auto trimmed_angular_velocity = findRoot(&trim, 2.0, 30.0);

	ac_input_state.rotor_inputs[rotor_idx].angular_velocity = trimmed_angular_velocity;

}

extern (C++) @nogc void no_trim(alias lift_model)(ref AircraftState ac_state, Aircraft aircraft, ref AircraftInputState ac_input_state, ref AircraftInputState last_ac_input_state, ref Atmosphere atmo, Inflow inflow, ref Wake wake, size_t rotor_idx, double time, double dt) {
	import std.conv : to;
	import std.math : PI;
	import std.numeric : findRoot;
	import std.stdio : writeln;

	immutable R = aircraft.rotors[rotor_idx].radius;
	immutable R_squared = R*R;
	immutable Omega = ac_input_state.rotor_inputs[rotor_idx].angular_velocity;

	//immutable desired_C_T = ac_input_state.rotor_inputs[rotor_idx].required_thrust/(atmo.density*PI*R_squared*R_squared*Omega*Omega);

	aircraft.rotors[rotor_idx].compute_rotor_properties!lift_model(
		ac_state.rotor_states[rotor_idx],
		ac_input_state.rotor_inputs[rotor_idx],
		ac_state,
		inflow,
		wake,
		ac_state.rotor_states[rotor_idx].C_T,
		time,
		dt
	);
}

string switchBuilder(int level, string switchVar, Args...)(string statement) {
	import std.conv : to;
	import std.string : indexOf, isNumeric;
	alias list = AliasSeq!Args;

	string fillPlaceHolder(int level)(in string statement, string arg)
	{
		auto strSlice = statement;
		ptrdiff_t searchIdx = 0;
		string newStatement = statement;
		while(searchIdx < statement.length)
		{
			auto idxStart = newStatement.indexOf('{', searchIdx);
			if(idxStart == -1)
			{
				break;
			}
			
			auto idxEnd = newStatement.indexOf('}', idxStart);
			assert(idxEnd > idxStart);
			
			auto sliceLen = idxEnd - (idxStart + 1);
			if(newStatement[idxStart+1..idxStart+1 + sliceLen].isNumeric)
			{
				auto strLevel = newStatement[idxStart+1..idxStart+1 + sliceLen].to!int;
				if(strLevel == level)
				{
					newStatement = newStatement[0..idxStart] ~ arg ~ newStatement[idxEnd + 1..$];
				}
			}

			searchIdx = idxStart + 1;
		}
		
		return newStatement;
	}

	string switchStatement = `final switch(`~switchVar~`)
	{`;
		foreach(arg; list)
		{
			auto thisStatement = fillPlaceHolder!level(statement, arg.stringof);
			switchStatement ~= `
	case `~arg.stringof~`:
		`~thisStatement~`
		/+break;+/`;
		}
		switchStatement ~= `
	}`;
	return switchStatement;
}

extern (C++) TrimFunction get_trim_with_lift(const char* _trim_algo, const char* _lift_model) {

	import std.conv : to;

	string trim_algo = _trim_algo.to!string;
	string lift_model = _lift_model.to!string;

	immutable str = `
			alias t = mixin({0});
			alias l = mixin({1});

			return &t!l;
		`
		.switchBuilder!(1, "lift_model", lift_model_list)
		.switchBuilder!(0, "trim_algo", trim_list)
	;

	mixin(str);
}
+/
