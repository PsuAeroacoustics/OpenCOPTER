module opencopter.aircraft.state;

import opencopter.aircraft;
import opencopter.config;
import opencopter.math;
import opencopter.memory;
import opencopter.weissingerl;
import opencopter.vortexlattice;

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
		immutable num_wings = ac.wings.length;
		immutable num_wing_parts = ac.wings[0].wing_parts.length;
		immutable num_wing_ctrl_chunks = ac.wings[0].wing_parts[0].ctrl_chunks.length;
		immutable num_span_chunks = ac.wings[0].wing_parts[0].chunks.length;
		immutable num_chord_nodes = num_wing_ctrl_chunks/num_span_chunks;
		mixin(array_ctor_mixin!(AC, "AircraftStateT!(AC)", "aircraft_history", "timesteps"));

		foreach(ref ac_hist; aircraft_history) {
			ac_hist = AircraftStateT!AC(num_rotors, num_blades, num_chunks*chunk_size, num_wings, num_wing_parts, num_span_chunks, num_chord_nodes, ac);
		}
	}

	this(AircraftT!AC* ac, size_t timesteps) {
		immutable num_blades = ac.rotors[0].blades.length;
		immutable num_rotors = ac.rotors.length;
		immutable num_chunks = ac.rotors[0].blades[0].chunks.length;
		immutable num_wings = ac.wings.length;
		immutable num_wing_parts = ac.wings[0].wing_parts.length;
		immutable num_wing_ctrl_chunks = ac.wings[0].wing_parts[0].ctrl_chunks.length;
		immutable num_span_chunks = ac.wings[0].wing_parts[0].chunks.length;
		immutable num_chord_nodes = num_wing_ctrl_chunks/num_span_chunks;
		mixin(array_ctor_mixin!(AC, "AircraftStateT!(AC)", "aircraft_history", "timesteps"));

		foreach(ref ac_hist; aircraft_history) {
			ac_hist = AircraftStateT!AC(num_rotors, num_blades, num_chunks*chunk_size, num_wings, num_wing_parts, num_span_chunks, num_chord_nodes, *ac);
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
	mixin ArrayDeclMixin!(AC, WingStateT!(AC), "wing_states");

	Vec4 freestream;

	this(size_t num_rotors, size_t num_blades, size_t num_elements, size_t num_wings, size_t num_wing_parts, size_t num_span_nodes, size_t num_chord_nodes, ref AircraftT!AC ac) {
		immutable actual_num_elements = num_elements%chunk_size == 0 ? num_elements : num_elements + (chunk_size - num_elements%chunk_size);

		immutable num_chunks = actual_num_elements/chunk_size;
		mixin(array_ctor_mixin!(AC, "RotorStateT!(AC)", "rotor_states", "num_rotors"));
		mixin(array_ctor_mixin!(AC, "WingStateT!(AC)", "wing_states", "num_wings"));
		foreach(i, ref rotor_state; rotor_states) {
			rotor_state = RotorStateT!AC(num_blades, num_chunks, ac.rotors[i]);
		}
		foreach(j, ref wing_state; wing_states) {
			wing_state = WingStateT!(AC)(num_wing_parts, num_span_nodes, num_chord_nodes, ac.wings[j]);
		}
	}

	this(size_t num_rotors, size_t num_blades, size_t num_elements, size_t num_wings, size_t num_wing_parts, size_t num_span_nodes, size_t num_chord_nodes, AircraftT!AC* ac) {
		immutable actual_num_elements = num_elements%chunk_size == 0 ? num_elements : num_elements + (chunk_size - num_elements%chunk_size);

		immutable num_chunks = actual_num_elements/chunk_size;
		mixin(array_ctor_mixin!(AC, "RotorStateT!(AC)", "rotor_states", "num_rotors"));
		mixin(array_ctor_mixin!(AC, "WingStateT!(AC)", "wing_states", "num_wings"));
		foreach(i, ref rotor_state; rotor_states) {
			rotor_state = RotorStateT!AC(num_blades, num_chunks, ac.rotors[i]);
		}
		foreach(j, ref wing_state; wing_states) {
			wing_state = WingStateT!(AC)(num_wing_parts, num_span_nodes, num_chord_nodes, ac.wings[j]);
		}
	}

	this(size_t num_rotors, size_t[] num_blades, size_t num_elements, size_t num_wings,size_t[] num_wing_parts,  size_t num_span_nodes, size_t num_chord_nodes, ref AircraftT!AC ac) {
		immutable actual_num_elements = num_elements%chunk_size == 0 ? num_elements : num_elements + (chunk_size - num_elements%chunk_size);

		immutable num_chunks = actual_num_elements/chunk_size;
		mixin(array_ctor_mixin!(AC, "RotorStateT!(AC)", "rotor_states", "num_rotors"));
		mixin(array_ctor_mixin!(AC, "WingStateT!(AC)", "wing_states", "num_wings"));
		foreach(i, ref rotor_state; rotor_states) {
			rotor_state = RotorStateT!AC(num_blades[i], num_chunks, ac.rotors[i]);
		}
		foreach(j, ref wing_state; wing_states) {
			wing_state = WingStateT!(AC)(num_wing_parts[j], num_span_nodes, num_chord_nodes, ac.wings[j]);
		}
	}

	this(size_t num_rotors, size_t[] num_blades, size_t num_elements, size_t num_wings,size_t[] num_wing_parts,  size_t num_span_nodes, size_t num_chord_nodes, AircraftT!AC* ac) {
		immutable actual_num_elements = num_elements%chunk_size == 0 ? num_elements : num_elements + (chunk_size - num_elements%chunk_size);

		immutable num_chunks = actual_num_elements/chunk_size;
		mixin(array_ctor_mixin!(AC, "RotorStateT!(AC)", "rotor_states", "num_rotors"));
		mixin(array_ctor_mixin!(AC, "WingStateT!(AC)", "wing_states", "num_wings"));
		foreach(i, ref rotor_state; rotor_states) {
			rotor_state = RotorStateT!AC(num_blades[i], num_chunks, ac.rotors[i]);
		}
		foreach(j, ref wing_state; wing_states) {
			wing_state = WingStateT!(AC)(num_wing_parts[j], num_span_nodes, num_chord_nodes, ac.wings[j]);
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

	mixin ArrayDeclMixin!(AC, WingPartStateT!(AC), "wing_part_states");
	double C_L;
	//double C_Mx;
	//double C_My;
	//double advance_ratio; // non-dim
	//double axial_advance_ratio; // non-dim

	this(size_t num_wing_parts, size_t span_elements, size_t chordwise_nodes, ref WingGeometryT!AC wing) {
		mixin(array_ctor_mixin!(AC, "WingPartStateT!(AC)", "wing_part_states", "num_wing_parts"));
		foreach(i, ref part_state; wing_part_states) {
			part_state = WingPartStateT!AC(span_elements, chordwise_nodes, wing.wing_parts[i]);
		}
	}

	this(size_t num_wing_parts, size_t span_elements, size_t chordwise_nodes, WingGeometryT!AC* wing) {
		mixin(array_ctor_mixin!(AC, "WingPartStateT!(AC)", "wing_part_states", "num_wing_parts"));
		foreach(i, ref part_state; wing_part_states) {
			part_state = WingPartStateT!AC(span_elements, chordwise_nodes, wing.wing_parts[i]);
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
		this.C_L = wing.C_L;
		return this;
	}

	ref typeof(this) opAssign(ref typeof(this) wing) {
		import std.stdio : writeln;
		//debug writeln("BladeGeometryT opAssign");
		this.wing_part_states = wing.wing_part_states;
		this.C_L = wing.C_L;
		return this;
	}

	ref typeof(this) opAssign(typeof(this)* wing) {
		import std.stdio : writeln;
		//debug writeln("BladeGeometryT opAssign");
		this.wing_part_states = wing.wing_part_states;
		this.C_L = wing.C_L;
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
	
	Chunk ctrl_pt_up;

	Chunk ctrl_pt_ut;
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
	mixin ArrayDeclMixin!(AC,WingPartCtrlPointStateChunk, "ctrl_chunks");

	VortexLatticeT!AC* circulation_model;

	this(size_t span_elements,size_t chord_ctrl_pt, ref WingPartGeometryT!AC wing_part) {
		immutable actual_num_span_elements = span_elements%chunk_size == 0 ? span_elements : span_elements + (chunk_size - span_elements%chunk_size);
		//enforce(num_elements % chunk_size == 0, "Number of spanwise elements must be a multiple of the chunk size ("~chunk_size.to!string~")");
		immutable span_chunks = actual_num_span_elements/chunk_size;
		size_t ctrl_pt_chunks = span_chunks*chord_ctrl_pt;
		mixin(array_ctor_mixin!(AC, "WingPartStateChunk", "chunks", "span_chunks"));
		mixin(array_ctor_mixin!(AC, "WingPartCtrlPointStateChunk", "ctrl_chunks", "ctrl_pt_chunks"));
		circulation_model = new VortexLatticeT!AC(span_chunks*chunk_size, chord_ctrl_pt, wing_part);
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

	this(size_t span_elements,size_t chord_ctrl_pt, WingPartGeometryT!AC* wing_part) {
		assert(wing_part !is null);
		immutable actual_num_span_elements = span_elements%chunk_size == 0 ? span_elements : span_elements + (chunk_size - span_elements%chunk_size);
		//enforce(num_elements % chunk_size == 0, "Number of spanwise elements must be a multiple of the chunk size ("~chunk_size.to!string~")");
		immutable span_chunks = actual_num_span_elements/chunk_size;
		size_t ctrl_pt_chunks = span_chunks*chord_ctrl_pt;
		mixin(array_ctor_mixin!(AC, "WingPartStateChunk", "chunks", "span_chunks"));
		mixin(array_ctor_mixin!(AC, "WingPartCtrlPointStateChunk", "ctrl_chunks", "ctrl_pt_chunks"));
		circulation_model = new VortexLatticeT!AC(span_chunks*chunk_size, chord_ctrl_pt, *wing_part);
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

	ref typeof(this) opAssign(ref typeof(this) wing_part) {
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

	ref typeof(this) opAssign(typeof(this)* wing_part) {
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

double[][] get_wing_state_matrix(string value, ArrayContainer AC)(ref WingPartStateT!AC wing_part){
	
	size_t spanwise_chunks = wing_part.chunks.length;
	size_t chordwise_nodes = wing_part.ctrl_chunks.length/spanwise_chunks;
	double[][] state_matrix = allocate_dense(chordwise_nodes, spanwise_chunks*chunk_size);

	foreach(c_idx, ref ctrl_chunk; wing_part.ctrl_chunks){
		immutable sp_idx = c_idx < spanwise_chunks ? c_idx : c_idx%spanwise_chunks;
		immutable chrd_idx = (c_idx - c_idx%spanwise_chunks)/spanwise_chunks;

		immutable out_start_idx = sp_idx*chunk_size;
		
		immutable remaining = spanwise_chunks*chunk_size - out_start_idx;

		immutable out_end_idx = remaining > chunk_size ? (sp_idx +1)*chunk_size : out_start_idx + remaining;
		immutable in_end_idx = remaining > chunk_size ? chunk_size : remaining;

		mixin("state_matrix[chrd_idx][out_start_idx..out_end_idx] = ctrl_chunk."~value~"[0..in_end_idx];");
	}
	return state_matrix;
}


double[] get_wing_state_array(string value, ArrayContainer AC)(ref WingStateT!AC wing){
	
	size_t elements = wing.wing_part_states[0].chunks.length*chunk_size;
	double[] state_array = new double[2*elements];

	foreach(wp_idx, wing_part; wing.wing_part_states){
		foreach(c_idx, chunk; wing_part.chunks){
			
			immutable out_start_idx = wp_idx*elements + c_idx*chunk_size;
			
			immutable remaining = elements- out_start_idx;

			immutable out_end_idx = remaining > chunk_size ? wp_idx*elements + (c_idx + 1)*chunk_size : wp_idx*elements + out_start_idx + remaining;
			immutable in_end_idx = remaining > chunk_size ? chunk_size : remaining;

			mixin("state_array[out_start_idx..out_end_idx] = chunk."~value~"[0..in_end_idx];");
		}
	}
	return state_array;
}
