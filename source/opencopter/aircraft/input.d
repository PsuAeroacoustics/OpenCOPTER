module opencopter.aircraft.input;

import opencopter.aircraft;
import opencopter.config;
import opencopter.memory;

import std.traits;

template is_aircraft_input_state(A) {
	enum bool is_aircraft_input_state = {
		static if(isPointer!(A)) {
			return isInstanceOf!(AircraftInputStateT, PointerTarget!A);
		} else {
			return isInstanceOf!(AircraftInputStateT, A);
		}
	}();
}

alias AircraftInputState = AircraftInputStateT!(ArrayContainer.none);

struct AircraftInputStateT(ArrayContainer AC) {
	mixin ArrayDeclMixin!(AC, RotorInputStateT!(AC), "rotor_inputs");

	this(size_t num_rotors, size_t num_blades) {
		mixin(array_ctor_mixin!(AC, "RotorInputStateT!(AC)", "rotor_inputs", "num_rotors"));
		
		mixin(array_ctor_mixin!(AC, "WingInputStateT!(AC)", "wing_inputs", "num_wings"))
		
		foreach(ref rotor; rotor_inputs) {
			mixin(array_ctor_mixin!(AC, "double", "rotor.r_0", "num_blades"));
			mixin(array_ctor_mixin!(AC, "double", "rotor.blade_pitches", "num_blades"));
			mixin(array_ctor_mixin!(AC, "double", "rotor.blade_flapping_rate", "num_blades"));
			mixin(array_ctor_mixin!(AC, "double", "rotor.blade_flapping", "num_blades"));
		}
	}

	this(size_t num_rotors, size_t[] num_blades) {
		mixin(array_ctor_mixin!(AC, "RotorInputStateT!(AC)", "rotor_inputs", "num_rotors"));
		
		mixin(array_ctor_mixin!(AC, "WingInputStateT!(AC)", "wing_inputs", "num_wings"))
		
		foreach(r_idx, ref rotor; rotor_inputs) {
			mixin(array_ctor_mixin!(AC, "double", "rotor.r_0", "num_blades[r_idx]"));
			mixin(array_ctor_mixin!(AC, "double", "rotor.blade_pitches", "num_blades[r_idx]"));
			mixin(array_ctor_mixin!(AC, "double", "rotor.blade_flapping_rate", "num_blades[r_idx]"));
			mixin(array_ctor_mixin!(AC, "double", "rotor.blade_flapping", "num_blades[r_idx]"));
		}
	}
	// Orientation data?
}

template is_rotor_input_state(A) {
	enum bool is_rotor_input_state = {
		static if(isPointer!(A)) {
			return isInstanceOf!(RotorInputStateT, PointerTarget!A);
		} else {
			return isInstanceOf!(RotorInputStateT, A);
		}
	}();
}

alias RotorInputState = RotorInputStateT!(ArrayContainer.none);
extern (C++) struct RotorInputStateT(ArrayContainer AC) {
	double angular_velocity; // rad/s
	double angular_accel; // rad/s^2
	double azimuth; // radians

	mixin ArrayDeclMixin!(AC, double, "r_0");
	mixin ArrayDeclMixin!(AC, double, "blade_pitches");
	mixin ArrayDeclMixin!(AC, double, "blade_flapping_rate");
	mixin ArrayDeclMixin!(AC, double, "blade_flapping");
}

template is_wing_input_state(A) {
	enum bool is_wing_input_state = {
		static if(isPointer!(A)) {
			return isInstanceOf!(WingInputStateT, PointerTarget!A);
		} else {
			return isInstanceOf!(WingInputStateT, A);
		}
	}();
}

alias WingInputState = WingInputStateT!(ArrayContainer.none);
extern (C++) struct WingInputStateT(ArrayContainer AC) {
	double angle_of_attack; // rad
	double freestream_velocity; // m/s
	//double angular_velocity; // rad/s
	//double angular_accel; // rad/s^2
	//double azimuth; // radians

	//mixin ArrayDeclMixin!(AC, double, "r_0");
	//mixin ArrayDeclMixin!(AC, double, "blade_pitches");
	//mixin ArrayDeclMixin!(AC, double, "blade_flapping_rate");
	//mixin ArrayDeclMixin!(AC, double, "blade_flapping");
}