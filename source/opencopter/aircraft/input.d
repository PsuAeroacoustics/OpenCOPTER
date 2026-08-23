module opencopter.aircraft.input;

import opencopter.aircraft;
import opencopter.config;
import opencopter.memory;

import std.math : abs, fmod, PI, sgn;
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
	mixin ArrayDeclMixin!(AC, WingInputStateT!(AC), "wing_inputs");

	this(size_t num_rotors, size_t num_blades, size_t num_wings, size_t num_chunks = 0) {
		mixin(array_ctor_mixin!(AC, "RotorInputStateT!(AC)", "rotor_inputs", "num_rotors"));
		
		mixin(array_ctor_mixin!(AC, "WingInputStateT!(AC)", "wing_inputs", "num_wings"));
		
		foreach(r_idx, ref rotor; rotor_inputs) {
			mixin(array_ctor_mixin!(AC, "double", "rotor.r_0", "num_blades"));
			mixin(array_ctor_mixin!(AC, "double", "rotor.blade_pitches", "num_blades"));
			mixin(array_ctor_mixin!(AC, "double", "rotor.blade_flapping_rate", "num_blades"));
			mixin(array_ctor_mixin!(AC, "double", "rotor.blade_flapping", "num_blades"));
			mixin(array_ctor_mixin!(AC, "BladeInputStateT!(AC)", "rotor.blade_inputs", "num_blades"));
			if(num_chunks > 0) {
				foreach(b_idx; 0..rotor.blade_inputs.length)
					rotor.blade_inputs[b_idx] = BladeInputStateT!(AC)(num_chunks);
			}
		}
	}

	this(size_t num_rotors, size_t[] num_blades, size_t num_wings) {
		mixin(array_ctor_mixin!(AC, "RotorInputStateT!(AC)", "rotor_inputs", "num_rotors"));
		
		mixin(array_ctor_mixin!(AC, "WingInputStateT!(AC)", "wing_inputs", "num_wings"));

		foreach(r_idx, ref rotor; rotor_inputs) {
			mixin(array_ctor_mixin!(AC, "double", "rotor.r_0", "num_blades[r_idx]"));
			mixin(array_ctor_mixin!(AC, "double", "rotor.blade_pitches", "num_blades[r_idx]"));
			mixin(array_ctor_mixin!(AC, "double", "rotor.blade_flapping_rate", "num_blades[r_idx]"));
			mixin(array_ctor_mixin!(AC, "double", "rotor.blade_flapping", "num_blades[r_idx]"));
			mixin(array_ctor_mixin!(AC, "BladeInputStateT!(AC)", "rotor.blade_inputs", "num_blades[r_idx]"));
		}
	}

	/++
	 +	Constructor that also allocates per-station deflection/velocity
	 +	arrays on each `BladeInputStateT`. `num_chunks[r_idx]` is the
	 +	number of radial stations for rotor `r_idx`. Pass 0 for any rotor
	 +	that should keep empty per-station arrays.
	 +/
	this(size_t num_rotors, size_t[] num_blades, size_t num_wings, size_t[] num_chunks) {
		this(num_rotors, num_blades, num_wings);
		foreach(r_idx, ref rotor; rotor_inputs) {
			if(num_chunks.length > r_idx && num_chunks[r_idx] > 0) {
				foreach(b_idx; 0..rotor.blade_inputs.length)
					rotor.blade_inputs[b_idx] = BladeInputStateT!(AC)(num_chunks[r_idx]);
			}
		}
	}
	// Orientation data?
}

template is_blade_input_state(A) {
	enum bool is_blade_input_state = {
		static if(isPointer!(A)) {
			return isInstanceOf!(BladeInputStateT, PointerTarget!A);
		} else {
			return isInstanceOf!(BladeInputStateT, A);
		}
	}();
}

alias BladeInputState = BladeInputStateT!(ArrayContainer.none);

/++
 +	Per-station blade deflection input state.
 +
 +	Contains per-station deflection arrays that represent deltas from
 +	the base blade shape. These deflections are combined with the base
 +	geometry when computing global x,y,z positions of each radial station.
 +
 +	Uses plain Chunk[] arrays (not ArrayDeclMixin) so the per-station
 +	data is always accessible via value semantics (chunk[idx][]), which
 +	is required for the SIMD-style position/twist accumulation in the
 +	blade element solver. This holds for both ArrayContainer.none and
 +	ArrayContainer.array variants.
 +/
extern (C++) struct BladeInputStateT(ArrayContainer AC) {
	/++
	 +	Per-blade collective pitch (radians).
	 +	Default: `inf` (sentinel meaning "not set", fall back to legacy array).
	 +/
	double pitch = double.infinity;
	/++
	 +	Per-blade flapping angle (radians).
	 +	Default: `inf` (sentinel meaning "not set", fall back to legacy array).
	 +/
	double flapping = double.infinity;
	/++
	 +	Per-blade flapping rate (rad/s).
	 +	Default: `inf` (sentinel meaning "not set", fall back to legacy array).
	 +/
	double flapping_rate = double.infinity;
	/++
	 +	Per-blade initial vortex core size (non-dim).
	 +	Default: `inf` (sentinel meaning "not set", fall back to legacy array).
	 +/
	double r_0 = double.infinity;

	/++
	 +	Per-station flapping displacement (non-dim, z direction in blade frame).
	 +/
	Chunk[] flap_deflection;
	/++
	 +	Per-station lagging displacement (non-dim, x direction in blade frame).
	 +/
	Chunk[] lag_deflection;
	/++
	 +	Per-station twist angle delta (radians).
	 +/
	Chunk[] twist_deflection;
	/++
	 +	Per-station flapping velocity (non-dim, z direction in blade frame).
	 +/
	Chunk[] flap_velocity;
	/++
	 +	Per-station lagging velocity (non-dim, x direction in blade frame).
	 +/
	Chunk[] lag_velocity;

	this(size_t num_chunks) {
		flap_deflection = new Chunk[num_chunks];
		lag_deflection = new Chunk[num_chunks];
		twist_deflection = new Chunk[num_chunks];
		flap_velocity = new Chunk[num_chunks];
		lag_velocity = new Chunk[num_chunks];
	}
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

	/++
	 +	DEPRECATED: Use `blade_inputs[blade_idx].pitch` instead.
	 +	Per-blade collective pitch (radians).
	 +/
	mixin ArrayDeclMixin!(AC, double, "r_0");
	/++
	 +	DEPRECATED: Use `blade_inputs` for new code.
	 +	Per-blade pitch (radians), read every timestep.
	 +/
	mixin ArrayDeclMixin!(AC, double, "blade_pitches");
	/++
	 +	DEPRECATED: Use `blade_inputs` for new code.
	 +	Per-blade flapping rate (rad/s), read every timestep.
	 +/
	mixin ArrayDeclMixin!(AC, double, "blade_flapping_rate");
	/++
	 +	DEPRECATED: Use `blade_inputs` for new code.
	 +	Per-blade flapping angle (radians), read every timestep.
	 +/
	mixin ArrayDeclMixin!(AC, double, "blade_flapping");

	/++
	 +	Per-blade deflection input states (new API).
	 +	Each element contains per-station flap/lag/twist deflection arrays.
	 +/
	mixin ArrayDeclMixin!(AC, BladeInputStateT!(AC), "blade_inputs");
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
	double cos_aoa;
	double sin_aoa;
	double freestream_velocity; // m/s
}

void basic_aircraft_rotor_dynamics(RIS)(auto ref RIS ac_input, double dt) {
	foreach(r_idx, ref rotor; ac_input.rotor_inputs) {
		rotor.azimuth += rotor.angular_velocity*dt + rotor.angular_accel*dt*dt;
		auto sign = sgn(rotor.azimuth);
		// Keep the azimuth between 0 and 2*PI so we don't
		// lose fp precicion as the sim marches in time and
		// the azimuth grows unbounded.
		if(abs(rotor.azimuth) > 2.0*PI) {
			rotor.azimuth = sign * fmod(abs(rotor.azimuth), 2.0*PI);
		}
	}
}

double basic_single_rotor_dynamics(RIS)(auto ref RIS input_state, double dt) {
	double angle = input_state.angular_velocity*dt + input_state.angular_accel*dt*dt;

	input_state.azimuth += angle;

	// Keep the azimuth between 0 and 2*PI so we don't
	// lose fp precicion as the sim marches in time and
	// the azimuth grows unbounded.
	if(input_state.azimuth > 2.0*PI) {
		input_state.azimuth = fmod(abs(input_state.azimuth), 2.0*PI);
	}

	return angle;
}
