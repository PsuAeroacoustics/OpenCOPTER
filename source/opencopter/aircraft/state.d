module opencopter.aircraft.state;

import opencopter.aircraft;
import opencopter.config;
import opencopter.memory;
import opencopter.weissingerl;

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

	this(size_t num_rotors, size_t num_blades, size_t num_elements, size_t num_wings, size_t num_wing_parts, size_t num_span_chunks, size_t num_chord_nodes, ref AircraftT!AC ac) {
		immutable actual_num_elements = num_elements%chunk_size == 0 ? num_elements : num_elements + (chunk_size - num_elements%chunk_size);
		//enforce(num_elements % chunk_size == 0, "Number of spanwise elements must be a multiple of the chunk size ("~chunk_size.to!string~")");
		immutable num_chunks = actual_num_elements/chunk_size;
		mixin(array_ctor_mixin!(AC, "RotorStateT!(AC)", "rotor_states", "num_rotors"));
		mixin(array_ctor_mixin!(AC, "WingStateT!(AC)", "wing_states", "num_wings"));
		foreach(i, ref rotor_state; rotor_states) {
			rotor_state = RotorStateT!AC(num_blades, num_chunks, ac.rotors[i]);
		}
		foreach(j, ref wing_state; wing_states) {
			wing_state = WingStateT!(AC)(num_wing_parts, num_span_chunks, num_chord_nodes, ac.wings[j]);
		}
	}

	this(size_t num_rotors, size_t num_blades, size_t num_elements, size_t num_wings, size_t num_wing_parts, size_t num_span_chunks, size_t num_chord_nodes, AircraftT!AC* ac) {
		immutable actual_num_elements = num_elements%chunk_size == 0 ? num_elements : num_elements + (chunk_size - num_elements%chunk_size);
		//enforce(num_elements % chunk_size == 0, "Number of spanwise elements must be a multiple of the chunk size ("~chunk_size.to!string~")");
		immutable num_chunks = actual_num_elements/chunk_size;
		mixin(array_ctor_mixin!(AC, "RotorStateT!(AC)", "rotor_states", "num_rotors"));
		mixin(array_ctor_mixin!(AC, "WingStateT!(AC)", "wing_states", "num_wings"));
		foreach(i, ref rotor_state; rotor_states) {
			rotor_state = RotorStateT!AC(num_blades, num_chunks, ac.rotors[i]);
		}
		foreach(j, ref wing_state; wing_states) {
			wing_state = WingStateT!(AC)(num_wing_parts, num_span_chunks, num_chord_nodes, ac.wings[j]);
		}
	}

	this(size_t num_rotors, size_t[] num_blades, size_t num_elements, size_t num_wings,size_t[] num_wing_parts,  size_t num_span_chunks, size_t num_chord_nodes, ref AircraftT!AC ac) {
		immutable actual_num_elements = num_elements%chunk_size == 0 ? num_elements : num_elements + (chunk_size - num_elements%chunk_size);
		//enforce(num_elements % chunk_size == 0, "Number of spanwise elements must be a multiple of the chunk size ("~chunk_size.to!string~")");
		immutable num_chunks = actual_num_elements/chunk_size;
		mixin(array_ctor_mixin!(AC, "RotorStateT!(AC)", "rotor_states", "num_rotors"));
		mixin(array_ctor_mixin!(AC, "WingStateT!(AC)", "wing_states", "num_wings"));
		foreach(i, ref rotor_state; rotor_states) {
			rotor_state = RotorStateT!AC(num_blades[i], num_chunks, ac.rotors[i]);
		}
		foreach(j, ref wing_state; wing_states) {
			wing_state = WingStateT!(AC)(num_wing_parts[j], num_span_chunks, num_chord_nodes, ac.wings[j]);
		}
	}

	this(size_t num_rotors, size_t[] num_blades, size_t num_elements, size_t num_wings,size_t[] num_wing_parts,  size_t num_span_chunks, size_t num_chord_nodes, AircraftT!AC* ac) {
		immutable actual_num_elements = num_elements%chunk_size == 0 ? num_elements : num_elements + (chunk_size - num_elements%chunk_size);
		//enforce(num_elements % chunk_size == 0, "Number of spanwise elements must be a multiple of the chunk size ("~chunk_size.to!string~")");
		immutable num_chunks = actual_num_elements/chunk_size;
		mixin(array_ctor_mixin!(AC, "RotorStateT!(AC)", "rotor_states", "num_rotors"));
		mixin(array_ctor_mixin!(AC, "WingStateT!(AC)", "wing_states", "num_wings"));
		foreach(i, ref rotor_state; rotor_states) {
			rotor_state = RotorStateT!AC(num_blades[i], num_chunks, ac.rotors[i]);
		}
		foreach(j, ref wing_state; wing_states) {
			wing_state = WingStateT!(AC)(num_wing_parts[j], num_span_chunks, num_chord_nodes, ac.wings[j]);
		}
	}

	ref typeof(this) opAssign(typeof(this) ac) {
		this.rotor_states = ac.rotor_states;
		this.wing_states = ac.wing_states;
		return this;
	}

	ref typeof(this) opAssign(ref typeof(this) ac) {
		this.rotor_states = ac.rotor_states;
		this.wing_states = ac.wing_states;
		return this;
	}

	ref typeof(this) opAssign(typeof(this)* ac) {
		this.rotor_states = ac.rotor_states;
		this.wing_states = ac.wing_states;
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
	double C_Mx;
	double C_My;
	double advance_ratio; // non-dim
	double axial_advance_ratio; // non-dim

	this(size_t num_blades, size_t num_chunks, ref RotorGeometryT!AC rotor) {
		mixin(array_ctor_mixin!(AC, "BladeStateT!(AC)", "blade_states", "num_blades"));
		foreach(i, ref blade_state; blade_states) {
			blade_state = BladeStateT!AC(num_chunks, rotor.blades[i]);
		}
	}

	this(size_t num_blades, size_t num_chunks, RotorGeometryT!AC* rotor) {
		mixin(array_ctor_mixin!(AC, "BladeStateT!(AC)", "blade_states", "num_blades"));
		foreach(i, ref blade_state; blade_states) {
			blade_state = BladeStateT!AC(num_chunks, rotor.blades[i]);
		}
	}

	@nogc ~this() {
		import std.stdio : writeln;
		//debug writeln("Destroying RotorState");
	}

	ref typeof(this) opAssign(typeof(this) rotor) {
		import std.stdio : writeln;
		//debug writeln("BladeGeometryT opAssign");
		this.blade_states = rotor.blade_states;
		this.C_T = rotor.C_T;
		this.advance_ratio = rotor.advance_ratio;
		this.axial_advance_ratio = rotor.axial_advance_ratio;

		return this;
	}

	ref typeof(this) opAssign(ref typeof(this) rotor) {
		//import std.stdio : writeln;
		//debug writeln("BladeGeometryT ref opAssign: ", blade.chunks.length);
		this.blade_states = rotor.blade_states;
		this.C_T = rotor.C_T;
		this.advance_ratio = rotor.advance_ratio;
		this.axial_advance_ratio = rotor.axial_advance_ratio;
		return this;
	}

	ref typeof(this) opAssign(typeof(this)* rotor) {
		//import std.stdio : writeln;
		//debug writeln("BladeGeometryT ptr opAssign: ", blade.chunks.length);
		this.blade_states = rotor.blade_states;
		this.C_T = rotor.C_T;
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
	 +  Spanwise lift coefficient distribution time derivative
	 +/
	Chunk dC_L_dot;
	/++
	 +  Spanwise drag coefficient distribution
	 +/
	Chunk dC_D;
	Chunk dC_Mx;
	Chunk dC_My;
	Chunk u_p;
	Chunk u_t;
	/++
	 +	Spanwise angle of attack
	 +/
	Chunk aoa;
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

	this(size_t num_chunks, ref BladeGeometryT!AC blade) {
		mixin(array_ctor_mixin!(AC, "BladeStateChunk", "chunks", "num_chunks"));

		circulation_model = new WeissingerL!AC(num_chunks*chunk_size, blade);
		foreach(ref chunk; chunks) {
			chunk.dC_L[] = 0;
			chunk.dC_T[] = 0;
			chunk.dC_D[] = 0;
			chunk.dC_Q[] = 0;
			chunk.aoa[] = 0;
			chunk.gamma[] = 0;
			chunk.d_gamma[] = 0;
			chunk.dC_Mx[] = 0;
			chunk.dC_My[] = 0;
		}
	}

	this(size_t num_chunks, BladeGeometryT!AC* blade) {
		assert(blade !is null);
		mixin(array_ctor_mixin!(AC, "BladeStateChunk", "chunks", "num_chunks"));

		circulation_model = new WeissingerL!AC(num_chunks*chunk_size, *blade);
		foreach(ref chunk; chunks) {
			chunk.dC_L[] = 0;
			chunk.dC_T[] = 0;
			chunk.dC_D[] = 0;
			chunk.dC_Q[] = 0;
			chunk.aoa[] = 0;
			chunk.gamma[] = 0;
			chunk.d_gamma[] = 0;
			chunk.dC_Mx[] = 0;
			chunk.dC_My[] = 0;
		}
	}

	ref typeof(this) opAssign(typeof(this) blade) {
		this.chunks = blade.chunks;
		this.azimuth = blade.azimuth;
		this.C_T = blade.C_T;
		this.C_Q = blade.C_Q;
		this.C_L = blade.C_L;
		this.C_D = blade.C_D;
		this.C_Mx = blade.C_Mx;
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
		this.C_Mx = blade.C_Mx;
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
		this.C_Mx = blade.C_Mx;
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
	double C_Mx;
	double C_My;
}


//Adding wing states 

template is_wing_state(A) {
	enum bool is_wing_state = {
		static if(isPointer!(A)) {
			return isInstanceOf!(WingStateT, PointerTarget!A);
		} else {
			return isInstanceOf!(WingStateT, A);
		}
	}();
}

alias WingState = WingStateT!(ArrayContainer.none);

extern (C++) struct WingStateT(ArrayContainer AC) {

	mixin ArrayDeclMixin!(AC, WingStateT!(AC), "wing_states");
	double C_T;
	//double C_Mx;
	//double C_My;
	//double advance_ratio; // non-dim
	//double axial_advance_ratio; // non-dim

	this(size_t num_wing_parts, size_t num_span_chunks, size_t chordwise_nodes, ref WingGeometryT!AC wing) {
		mixin(array_ctor_mixin!(AC, "WingPartStateT!(AC)", "wing_part_states", "num_parts"));
		foreach(i, ref part_state; wing_part_states) {
			wing_state = WingPartStateT!AC(num_span_chunks, chordwise_nodes, wing.wing_parts[i]);
		}
	}

	this(size_t num_wing_parts, size_t num_span_chunks, size_t chordwise_nodes, ref WingGeometryT!AC* wing) {
		mixin(array_ctor_mixin!(AC, "WingPartStateT!(AC)", "wing_part_states", "num_parts"));
		foreach(i, ref part_state; wing_part_states) {
			wing_part_state = WingPartStateT!AC(num_span_chunks, chordwise_nodes, wing.wing_parts[i]);
		}
	}

	@nogc ~this() {
		import std.stdio : writeln;
		//debug writeln("Destroying RotorState");
	}

	ref typeof(this) opAssign(typeof(this) wing) {
		import std.stdio : writeln;
		//debug writeln("BladeGeometryT opAssign");
		this.wing_part_states = wing.wing_part_states;
		this.C_T = wing.C_T;
		return this;
	}

	ref typeof(this) opAssign(ref typeof(this) wing) {
		import std.stdio : writeln;
		//debug writeln("BladeGeometryT opAssign");
		this.wing_part_states = wing.wing_part_states;
		this.C_T = wing.C_T;
		return this;
	}

	ref typeof(this) opAssign(typeof(this)* wing) {
		import std.stdio : writeln;
		//debug writeln("BladeGeometryT opAssign");
		this.wing_part_states = wing.wing_part_states;
		this.C_T = wing.C_T;
		return this;
	}
}


template is_wing_part_state_chunk(A) {
	enum bool is_Wing_Part_state_chunk = {
		static if(isPointer!(A)) {
			return isInstanceOf!(WingPartStateChunk, PointerTarget!A);
		} else {
			return isInstanceOf!(WingPartStateChunk, A);
		}
	}();
}

extern (C++) struct WingPartStateChunk{
	
	// Distribution of circulation
	//Chunk gamma;

	// Distribution of spanwise change in circulation 
	//Chunk d_gamma;

	// distribution of elemental angle of attack
	Chunk aoa;

	// distribution of coefficient of lift
	Chunk dC_L;

	// distribution of coefficient of pitching moment
	Chunk dC_M;
}

template is_wing_part_ctrl_pt_state_chunk(A) {
	enum bool is_Wing_Part_ctrl_pt_state_chunk = {
		static if(isPointer!(A)) {
			return isInstanceOf!(WingPartCtrlPointStateChunk, PointerTarget!A);
		} else {
			return isInstanceOf!(WingPartCtrlPointStateChunk, A);
		}
	}();
}

extern (C++) struct WingPartCtrlPointStateChunk{
	// Distribution of angle of attack at each control point
	Chunk ctrl_pt_aoa;
}

/*extern (C++) struct WingCirculationChunk{
	// Distribution of circulation
	Chunk gamma;

	// Distribution of spanwise change in circulation 
	//Chunk d_gamma;
}*/

template is_wing_part_state(A) {
	enum bool is_wing_part_state = {
		static if(isPointer!(A)) {
			return isInstanceOf!(WingPartStateT, PointerTarget!A);
		} else {
			return isInstanceOf!(WingPartStateT, A);
		}
	}();
}

alias WingPartState = WingPartStateT!(ArrayContainer.none);
/++
 +	Current state of the blade physics properties
 +/
extern (C++) struct WingPartStateT(ArrayContainer AC) {
	/++
	 +	Span wise blade data. Chunked together for cache locality and vectorization.
	 +/
	mixin ArrayDeclMixin!(AC,WingPartStateChunk, "chunks");

	VortexLattice!AC* circulation_model;

	this(size_t span_chunks,size_t chord_ctrl_pt, ref WingPartGeometryT!AC wing_part) {
		size_t ctrl_pt_chunks = span_chunks*chord_ctrl_pt;
		mixin(array_ctor_mixin!(AC, "WingPartStateChunk", "chunks", "span_chunks"));
		mixin(array_ctor_mixin!(AC, "WingPartCtrlPointStateChunk", "ctrl_chunks", "ctrl_pt_chunks"))
		circulation_model = new VortexLattice!AC(span_chunks*chunk_size, chord_ctrl_pt, wing_part);
		foreach(ref chunk; chunks) {
			chunk.dC_L[] = 0;
			//chunk.dC_T[] = 0;
			//chunk.dC_D[] = 0;
			//chunk.dC_Q[] = 0;
			chunk.aoa[] = 0;
			//chunk.gamma[] = 0;
			//chunk.d_gamma[] = 0;
			chunk.dC_M[] = 0;
			//chunk.dC_My[] = 0;
		}
		foreach(ref ctrl_chunk;ctrl_chunks){
			ctrl_chunk.ctrl_pt_aoa[]=0;
		}
	}

	this(size_t span_chunks,size_t chord_ctrl_pt, ref WingPartGeometryT!AC* wing_part) {
		assert(wing !is null);
		size_t ctrl_pt_chunks = span_chunks*chord_ctrl_pt;
		mixin(array_ctor_mixin!(AC, "WingPartStateChunk", "chunks", "span_chunks"));
		mixin(array_ctor_mixin!(AC, "WingPartCtrlPointStateChunk", "ctrl_chunks", "ctrl_pt_chunks"))
		circulation_model = new VortexLattice!AC(span_chunks*chunk_size, chord_ctrl_pt, *wing_part);
		foreach(ref chunk; chunks) {
			chunk.dC_L[] = 0;
			//chunk.dC_T[] = 0;
			//chunk.dC_D[] = 0;
			//chunk.dC_Q[] = 0;
			chunk.aoa[] = 0;
			//chunk.gamma[] = 0;
			//chunk.d_gamma[] = 0;
			chunk.dC_M[] = 0;
			//chunk.dC_My[] = 0;
		}
		foreach(ref ctrl_chunk;ctrl_chunks){
			ctrl_chunk.ctrl_pt_aoa[]=0;
		}
	}

	ref typeof(this) opAssign(typeof(this) wing_part) {
		this.chunks = wing_part.chunks;
		this.ctrl_chunks = wing_part.ctrl_chunks;
		//this.azimuth = blade.azimuth;
		//this.C_T = blade.C_T;
		//this.C_Q = blade.C_Q;
		this.C_L = wing_part.C_L;
		//this.C_D = blade.C_D;
		this.C_M = wing_part.C_M;
		//this.C_My = blade.C_My;
		this.circulation_model = wing_part.circulation_model;
		return this;
	}

	ref typeof(this) opAssign(ref typeof(this) wing) {
		this.chunks = wing_part.chunks;
		this.ctrl_chunks = wing_part.ctrl_chunks;
		//this.azimuth = blade.azimuth;
		//this.C_T = blade.C_T;
		//this.C_Q = blade.C_Q;
		this.C_L = wing_part.C_L;
		//this.C_D = blade.C_D;
		this.C_M = wing_part.C_M;
		//this.C_My = blade.C_My;
		this.circulation_model = wing_part.circulation_model;
		return this;
	}

	ref typeof(this) opAssign(typeof(this)* wing) {
		this.chunks = wing_part.chunks;
		this.ctrl_chunks = wing_part.ctrl_chunks;
		//this.azimuth = blade.azimuth;
		//this.C_T = blade.C_T;
		//this.C_Q = blade.C_Q;
		this.C_L = wing_part.C_L;
		//this.C_D = blade.C_D;
		this.C_M = wing_part.C_M;
		//this.C_My = blade.C_My;
		this.circulation_model = wing_part.circulation_model;
		return this;
	}

	//double azimuth;

	/++
	 +	Thrust coefficient for single blade.
	 +/
	//double C_T;
	/++
	 +  Power coefficient for single blade.
	 +/
	//double C_Q;
	/++
	 +  Lift coefficient for single blade.
	 +/
	double C_L;
	/++
	 +  Drag coefficient for single blade.
	 +/
	//double C_D;
	//double C_Mx;
	double C_M;
}

double[] get_state_array(string value, ArrayContainer AC)(ref WingStateT!AC wing) {
	immutable elements = wing.chunks.length*chunk_size;
	double[] state_array = new double[elements];
	foreach(c_idx, ref chunk; wing.chunks) {
		immutable out_start_idx = c_idx*chunk_size;

		immutable remaining = elements - out_start_idx;
		
		immutable out_end_idx = remaining > chunk_size ? (c_idx + 1)*chunk_size : out_start_idx + remaining;
		immutable in_end_idx = remaining > chunk_size ? chunk_size : remaining;

		mixin("state_array[out_start_idx..out_end_idx] = chunk."~value~"[0..in_end_idx];");
	}
	return state_array;
}

void get_state_array(string value, ArrayContainer AC)(ref WingStateT!AC wing, auto ref double[] state_array) {
	foreach(c_idx, ref chunk; wing.chunks) {
		immutable out_start_idx = c_idx*chunk_size;

		immutable remaining = state_array.length - out_start_idx;
		
		immutable out_end_idx = remaining > chunk_size ? (c_idx + 1)*chunk_size : out_start_idx + remaining;
		immutable in_end_idx = remaining > chunk_size ? chunk_size : remaining;

		mixin("state_array[out_start_idx..out_end_idx] = chunk."~value~"[0..in_end_idx];");
	}
}

void get_state_array(string value, ArrayContainer AC)(ref WingStateT!AC wing, auto ref float[] state_array) {
	foreach(c_idx, ref chunk; wing.chunks) {
		immutable out_start_idx = c_idx*chunk_size;

		immutable remaining = state_array.length - out_start_idx;
		
		immutable out_end_idx = remaining > chunk_size ? (c_idx + 1)*chunk_size : out_start_idx + remaining;
		immutable in_end_idx = remaining > chunk_size ? chunk_size : remaining;

		mixin("state_array[out_start_idx..out_end_idx] = chunk."~value~"[].map!(a => a.to!float).staticArray!(float[chunk_size])[0..in_end_idx];");
	}
}

void set_state_array(string value, ArrayContainer AC)(ref WingStateT!AC wing, double[] data) {

	foreach(c_idx, ref chunk; wing.chunks) {
		immutable out_start_idx = c_idx*chunk_size;

		immutable remaining = data.length - out_start_idx;
		
		immutable out_end_idx = remaining > chunk_size ? (c_idx + 1)*chunk_size : out_start_idx + remaining;
		immutable in_end_idx = remaining > chunk_size ? chunk_size : remaining;

		mixin("chunk."~value~"[0..in_end_idx] = data[out_start_idx..out_end_idx];");
	}
}
