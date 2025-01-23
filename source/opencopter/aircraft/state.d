module opencopter.aircraft.state;

import opencopter.aircraft;
import opencopter.config;
import opencopter.math;
import opencopter.memory;
import opencopter.weissingerl;

import numd.linearalgebra.matrix;

import std.algorithm : map;
import std.array : staticArray;
import std.conv : to;
import std.exception : enforce;
import std.traits;

alias AircraftTimehistory = AircraftTimehistoryT!(ArrayContainer.none);
extern (C++) struct AircraftTimehistoryT(ArrayContainer AC) {

	mixin ArrayDeclMixin!(AC, AircraftStateT!(AC), "aircraft_history");

	this(ref AircraftT!AC ac, size_t timesteps) {
		immutable num_blades = ac.rotors[0].blades.length;
		immutable num_rotors = ac.rotors.length;
		immutable num_chunks = ac.rotors[0].blades[0].chunks.length;
		mixin(array_ctor_mixin!(AC, "AircraftStateT!(AC)", "aircraft_history", "timesteps"));

		foreach(ref ac_hist; aircraft_history) {
			ac_hist = AircraftStateT!AC(num_rotors, num_blades, num_chunks*chunk_size, ac);
		}
	}

	this(AircraftT!AC* ac, size_t timesteps) {
		immutable num_blades = ac.rotors[0].blades.length;
		immutable num_rotors = ac.rotors.length;
		immutable num_chunks = ac.rotors[0].blades[0].chunks.length;
		mixin(array_ctor_mixin!(AC, "AircraftStateT!(AC)", "aircraft_history", "timesteps"));

		foreach(ref ac_hist; aircraft_history) {
			ac_hist = AircraftStateT!AC(num_rotors, num_blades, num_chunks*chunk_size, *ac);
		}
	}
}

template is_aircraft_state(A) {
	enum bool is_aircraft_state = {
		static if(isPointer!(A)) {
			return isInstanceOf!(AircraftStateT, PointerTarget!A);
		} else {
			return isInstanceOf!(AircraftStateT, A);
		}
	}();
}

alias AircraftState = AircraftStateT!(ArrayContainer.none);

 struct AircraftStateT(ArrayContainer _AC) {
	alias AC = _AC;
	mixin ArrayDeclMixin!(AC, RotorStateT!(AC), "rotor_states");

	Vec4 freestream;

	this(size_t num_rotors, size_t num_blades, size_t num_elements, ref AircraftT!AC ac) {
		immutable actual_num_elements = num_elements%chunk_size == 0 ? num_elements : num_elements + (chunk_size - num_elements%chunk_size);

		immutable num_chunks = actual_num_elements/chunk_size;
		mixin(array_ctor_mixin!(AC, "RotorStateT!(AC)", "rotor_states", "num_rotors"));
		foreach(i, ref rotor_state; rotor_states) {
			rotor_state = RotorStateT!AC(num_blades, num_chunks, ac.rotors[i]);
		}
	}

	this(size_t num_rotors, size_t num_blades, size_t num_elements, AircraftT!AC* ac) {
		immutable actual_num_elements = num_elements%chunk_size == 0 ? num_elements : num_elements + (chunk_size - num_elements%chunk_size);

		immutable num_chunks = actual_num_elements/chunk_size;
		mixin(array_ctor_mixin!(AC, "RotorStateT!(AC)", "rotor_states", "num_rotors"));
		foreach(i, ref rotor_state; rotor_states) {
			rotor_state = RotorStateT!AC(num_blades, num_chunks, ac.rotors[i]);
		}
	}

	this(size_t num_rotors, size_t[] num_blades, size_t num_elements, ref AircraftT!AC ac) {
		immutable actual_num_elements = num_elements%chunk_size == 0 ? num_elements : num_elements + (chunk_size - num_elements%chunk_size);

		immutable num_chunks = actual_num_elements/chunk_size;
		mixin(array_ctor_mixin!(AC, "RotorStateT!(AC)", "rotor_states", "num_rotors"));
		foreach(i, ref rotor_state; rotor_states) {
			rotor_state = RotorStateT!AC(num_blades[i], num_chunks, ac.rotors[i]);
		}
	}

	this(size_t num_rotors, size_t[] num_blades, size_t num_elements, AircraftT!AC* ac) {
		immutable actual_num_elements = num_elements%chunk_size == 0 ? num_elements : num_elements + (chunk_size - num_elements%chunk_size);

		immutable num_chunks = actual_num_elements/chunk_size;
		mixin(array_ctor_mixin!(AC, "RotorStateT!(AC)", "rotor_states", "num_rotors"));
		foreach(i, ref rotor_state; rotor_states) {
			rotor_state = RotorStateT!AC(num_blades[i], num_chunks, ac.rotors[i]);
		}
	}

	ref typeof(this) opAssign(typeof(this) ac) {
		this.rotor_states = ac.rotor_states;

		return this;
	}

	ref typeof(this) opAssign(ref typeof(this) ac) {
		this.rotor_states = ac.rotor_states;
		return this;
	}

	ref typeof(this) opAssign(typeof(this)* ac) {
		this.rotor_states = ac.rotor_states;
		return this;
	}
}

template is_rotor_state(A) {
	enum bool is_rotor_state = {
		static if(isPointer!(A)) {
			return isInstanceOf!(RotorStateT, PointerTarget!A);
		} else {
			return isInstanceOf!(RotorStateT, A);
		}
	}();
}

alias RotorState = RotorStateT!(ArrayContainer.none);

extern (C++) struct RotorStateT(ArrayContainer AC) {

	mixin ArrayDeclMixin!(AC, BladeStateT!(AC), "blade_states");
	double C_T;
	double C_Q;
	double C_Mx;
	double C_My;
	
	double advance_ratio; // non-dim
	double axial_advance_ratio; // non-dim

	this(size_t num_blades, size_t num_chunks, ref RotorGeometryT!AC rotor) {
		mixin(array_ctor_mixin!(AC, "BladeStateT!(AC)", "blade_states", "num_blades"));
		foreach(i, ref blade_state; blade_states) {
			blade_state = BladeStateT!AC(num_chunks, rotor.blades[i], rotor.radius);
		}
	}

	this(size_t num_blades, size_t num_chunks, RotorGeometryT!AC* rotor) {
		mixin(array_ctor_mixin!(AC, "BladeStateT!(AC)", "blade_states", "num_blades"));
		foreach(i, ref blade_state; blade_states) {
			blade_state = BladeStateT!AC(num_chunks, rotor.blades[i], rotor.radius);
		}
	}

	@nogc ~this() {
	}

	ref typeof(this) opAssign(typeof(this) rotor) {
		this.blade_states = rotor.blade_states;
		this.C_T = rotor.C_T;
		this.C_Q = rotor.C_Q;
		this.advance_ratio = rotor.advance_ratio;
		this.axial_advance_ratio = rotor.axial_advance_ratio;

		return this;
	}

	ref typeof(this) opAssign(ref typeof(this) rotor) {
		this.blade_states = rotor.blade_states;
		this.C_T = rotor.C_T;
		this.C_Q = rotor.C_Q;
		this.advance_ratio = rotor.advance_ratio;
		this.axial_advance_ratio = rotor.axial_advance_ratio;
		return this;
	}

	ref typeof(this) opAssign(typeof(this)* rotor) {
		this.blade_states = rotor.blade_states;
		this.C_T = rotor.C_T;
		this.C_Q = rotor.C_Q;
		this.advance_ratio = rotor.advance_ratio;
		this.axial_advance_ratio = rotor.axial_advance_ratio;
		return this;
	}
}

template is_blade_state_chunk(A) {
	enum bool is_blade_state_chunk = {
		static if(isPointer!(A)) {
			return is(PointerTarget!A == BladeStateChunk);
		} else {
			return is(A == BladeStateChunk);
		}
	}();
}

extern (C++) struct BladeStateChunk {
	/++
	 +  Spanwise thrust coefficient distribution
	 +/
	Chunk dC_T;
	/++
	 +  Spanwise drag coefficient distribution in blade frame
	 +/
	Chunk dC_Db;
	/++
	 +  Spanwise sectional normal force coefficient distribution
	 +/
	Chunk dC_N;
	/++
	 +  Spanwise sectional chordwise force coefficient distribution
	 +/
	Chunk dC_c;
	/++
	 +  Spanwise thrust coefficient distribution time derivative
	 +/
	Chunk dC_T_dot;
	/++
	 +  Spanwise torque coefficient distribution
	 +/
	Chunk dC_Q;
	/++
	 +  Spanwise lift coefficient distribution
	 +/
	Chunk dC_L;
	/++
	 +  Spanwise lift airfoil coefficient distribution
	 +/
	Chunk dC_l;
	/++
	 +  Spanwise lift coefficient distribution time derivative
	 +/
	Chunk dC_L_dot;
	/++
	 +  Spanwise drag coefficient distribution in AF frame
	 +/
	Chunk dC_D;
	/++
	 +  Spanwise airfoil drag coefficient distribution in AF frame
	 +/
	Chunk dC_d;
	
	Chunk dC_Mz;
	Chunk dC_My;
	//Chunk dC_M;
	Chunk u_p;
	Chunk shed_u_p;
	Chunk u_t;
	/++
	 +	Spanwise angle of attack
	 +/
	Chunk aoa;
	/++
	 +	Spanwise angle of attack
	 +/
	Chunk aoa_eff;
	/++
	 +	Spanwise angle of attack
	 +/
	Chunk inflow_angle;
	/++
	 +	Spanwise circulation
	 +/
	Chunk gamma;
	/++
	 +	Spanwise change in circulation
	 +/
	Chunk d_gamma;

	Chunk x;
	Chunk y;
	Chunk z;

	Chunk r_c;

	Chunk theta;

	Vector!(4, Chunk) projected_vel;

	Vector!(4, Chunk) blade_local_vel;
}

template is_blade_state(A) {
	enum bool is_blade_state = {
		static if(isPointer!(A)) {
			return isInstanceOf!(BladeStateT, PointerTarget!A);
		} else {
			return isInstanceOf!(BladeStateT, A);
		}
	}();
}

alias BladeState = BladeStateT!(ArrayContainer.none);
/++
 +	Current state of the blade physics properties
 +/
extern (C++) struct BladeStateT(ArrayContainer AC) {
	/++
	 +	Span wise blade data. Chunked together for cache locality and vectorization.
	 +/
	mixin ArrayDeclMixin!(AC, BladeStateChunk, "chunks");

	WeissingerL!AC* circulation_model;
	
	this(size_t num_chunks, ref BladeGeometryT!AC blade, double radius) {
		mixin(array_ctor_mixin!(AC, "BladeStateChunk", "chunks", "num_chunks"));

		circulation_model = new WeissingerL!AC(num_chunks*chunk_size, blade, radius);
		foreach(ref chunk; chunks) {
			chunk.aoa_eff[] = 0;
			chunk.dC_L[] = 0;
			chunk.dC_T[] = 0;
			chunk.dC_D[] = 0;
			chunk.dC_Q[] = 0;
			chunk.aoa[] = 0;
			chunk.gamma[] = 0;
			chunk.d_gamma[] = 0;
			chunk.dC_Mz[] = 0;
			chunk.dC_My[] = 0;
			chunk.r_c[] = 0;
		}
	}

	this(size_t num_chunks, BladeGeometryT!AC* blade, double radius) {
		assert(blade !is null);
		mixin(array_ctor_mixin!(AC, "BladeStateChunk", "chunks", "num_chunks"));

		circulation_model = new WeissingerL!AC(num_chunks*chunk_size, *blade, radius);
		foreach(ref chunk; chunks) {
			chunk.aoa_eff[] = 0;
			chunk.dC_L[] = 0;
			chunk.dC_T[] = 0;
			chunk.dC_D[] = 0;
			chunk.dC_Q[] = 0;
			chunk.aoa[] = 0;
			chunk.gamma[] = 0;
			chunk.d_gamma[] = 0;
			chunk.dC_Mz[] = 0;
			chunk.dC_My[] = 0;
			chunk.r_c[] = 0;
		}
	}

	ref typeof(this) opAssign(typeof(this) blade) {
		this.chunks = blade.chunks;
		this.azimuth = blade.azimuth;
		this.C_T = blade.C_T;
		this.C_Q = blade.C_Q;
		this.C_L = blade.C_L;
		this.C_D = blade.C_D;
		this.C_Mz = blade.C_Mz;
		this.C_My = blade.C_My;
		this.circulation_model = blade.circulation_model;
		return this;
	}

	ref typeof(this) opAssign(ref typeof(this) blade) {
		this.chunks = blade.chunks;
		this.azimuth = blade.azimuth;
		this.C_T = blade.C_T;
		this.C_Q = blade.C_Q;
		this.C_L = blade.C_L;
		this.C_D = blade.C_D;
		this.C_Mz = blade.C_Mz;
		this.C_My = blade.C_My;
		this.circulation_model = blade.circulation_model;
		return this;
	}

	ref typeof(this) opAssign(typeof(this)* blade) {
		this.chunks = blade.chunks;
		this.azimuth = blade.azimuth;
		this.C_T = blade.C_T;
		this.C_Q = blade.C_Q;
		this.C_L = blade.C_L;
		this.C_D = blade.C_D;
		this.C_Mz = blade.C_Mz;
		this.C_My = blade.C_My;
		this.circulation_model = blade.circulation_model;
		return this;
	}

	double azimuth;

	/++
	 +	Thrust coefficient for single blade.
	 +/
	double C_T;
	/++
	 +  Power coefficient for single blade.
	 +/
	double C_Q;
	/++
	 +  Lift coefficient for single blade.
	 +/
	double C_L;
	/++
	 +  Drag coefficient for single blade.
	 +/
	double C_D;
	double C_Mz;
	double C_My;
	//double C_M;
}

double[] get_state_array(string value, ArrayContainer AC)(ref BladeStateT!AC blade) {
	immutable elements = blade.chunks.length*chunk_size;
	double[] state_array = new double[elements];
	foreach(c_idx, ref chunk; blade.chunks) {
		immutable out_start_idx = c_idx*chunk_size;

		immutable remaining = elements - out_start_idx;
		
		immutable out_end_idx = remaining > chunk_size ? (c_idx + 1)*chunk_size : out_start_idx + remaining;
		immutable in_end_idx = remaining > chunk_size ? chunk_size : remaining;

		mixin("state_array[out_start_idx..out_end_idx] = chunk."~value~"[0..in_end_idx];");
	}
	return state_array;
}

void get_state_array(string value, ArrayContainer AC)(ref BladeStateT!AC blade, auto ref double[] state_array) {
	foreach(c_idx, ref chunk; blade.chunks) {
		immutable out_start_idx = c_idx*chunk_size;

		immutable remaining = state_array.length - out_start_idx;
		
		immutable out_end_idx = remaining > chunk_size ? (c_idx + 1)*chunk_size : out_start_idx + remaining;
		immutable in_end_idx = remaining > chunk_size ? chunk_size : remaining;

		mixin("state_array[out_start_idx..out_end_idx] = chunk."~value~"[0..in_end_idx];");
	}
}

void get_state_array(string value, ArrayContainer AC)(ref BladeStateT!AC blade, auto ref float[] state_array) {
	foreach(c_idx, ref chunk; blade.chunks) {
		immutable out_start_idx = c_idx*chunk_size;

		immutable remaining = state_array.length - out_start_idx;
		
		immutable out_end_idx = remaining > chunk_size ? (c_idx + 1)*chunk_size : out_start_idx + remaining;
		immutable in_end_idx = remaining > chunk_size ? chunk_size : remaining;

		mixin("state_array[out_start_idx..out_end_idx] = chunk."~value~"[].map!(a => a.to!float).staticArray!(float[chunk_size])[0..in_end_idx];");
	}
}

void set_state_array(string value, ArrayContainer AC)(ref BladeStateT!AC blade, double[] data) {

	foreach(c_idx, ref chunk; blade.chunks) {
		immutable out_start_idx = c_idx*chunk_size;

		immutable remaining = data.length - out_start_idx;
		
		immutable out_end_idx = remaining > chunk_size ? (c_idx + 1)*chunk_size : out_start_idx + remaining;
		immutable in_end_idx = remaining > chunk_size ? chunk_size : remaining;

		mixin("chunk."~value~"[0..in_end_idx] = data[out_start_idx..out_end_idx];");
	}
}
