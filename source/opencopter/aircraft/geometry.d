module opencopter.aircraft.geometry;

import opencopter.aircraft;
import opencopter.config;
import opencopter.math;
import opencopter.memory;
import opencopter.airfoilmodels;

import numd.linearalgebra.matrix;

import std.conv : to;
import std.exception : enforce;
import std.traits;
import std.typecons;
import std.stdio;

/+extern(C++) struct RadiusPoints {
	double[] r;
	size_t num_points;
}+/

double[] generate_radius_points(size_t n_sections) {
	import std.algorithm : map;
	import std.array : array;
	import std.math : cos, PI;
	import std.range : iota, retro;
	//RadiusPoints r;
	immutable num_points = n_sections%chunk_size == 0 ? n_sections : n_sections + (chunk_size - n_sections%chunk_size);
    return iota(1.0, num_points + 1.0).map!((n) {
    	immutable psi = n*PI/(num_points.to!double + 1.0);
    	auto r = 0.5*(cos(psi) + 1.0).to!double;
    	return r;
    }).retro.array;
}

/*double[] generate_spanwise_vortex_nodes(size_t span_n_sections) {
	import std.algorithm : map;
	import std.array : array;
	import std.math : cos, PI;
	import std.range : iota, retro;
	// Spanwise Votex nodes
	immutable num_points = span_n_sections%chunk_size == 0 ? span_n_sections : span_n_sections + (chunk_size - span_n_sections%chunk_size);
	return iota(1.0,num_points + 2.0).map!((l){
		immutable phi = (2*l-1)*PI/(2.0*(num_points.to!double + 1.0));
		auto span_vortex_pt = 0.25*(1-cos(phi)).to!double;
		return span_vortex_pt;
	}).array;
}*/


double[] generate_spanwise_control_points(size_t span_n_sections) {
	import std.algorithm : map;
	import std.array : array;
	import std.math : cos, PI;
	import std.range : iota, retro;
	// Spanwise Votex nodes
	immutable num_points = span_n_sections%chunk_size == 0 ? span_n_sections : span_n_sections + (chunk_size - span_n_sections%chunk_size);
	return iota(0.0,num_points).map!((j){
		immutable phi = j*PI/(num_points.to!double);
		auto span_ctrl_pt = 0.25*(1-cos(phi)).to!double;
		return span_ctrl_pt;
	}).array;
}

void set_wing_ctrl_pt_geometry(WG)(auto ref WG wing_geometry, size_t _spanwise_nodes, size_t _chordwise_nodes, double _camber){
    
    import std.math : cos,tan, PI;
    import std.math : abs;
	double wing_span = wing_geometry.wing_parts[0].wing_span; // half wing span
	Chunk span = 2*wing_span; // full wing span
	double root_chord = wing_geometry.wing_parts[0].wing_root_chord;
	size_t acutal_span_nodes = _spanwise_nodes%chunk_size == 0 ? _spanwise_nodes : _spanwise_nodes + (chunk_size - _spanwise_nodes%chunk_size);
	size_t _spanwise_chunks = acutal_span_nodes/chunk_size;
    auto span_ctrl_pt = generate_spanwise_control_points(acutal_span_nodes);
	auto chord_ctrl_pt = generate_chordwise_control_points(_chordwise_nodes);
    //auto span_ctrl_pt_left = generate_spanwise_control_points(_spanwise_chunks*chunk_size);
    //span_vr_nodes_left = span_vr_nodes_left.retro;
    //auto chord_vr_nodes = generate_chordwise_votex_nodes(_chordwise_nodes);

	if(wing_geometry.wing_parts.length == 1){
		string side = wing_geometry.wing_parts[0].loc;
		writeln("location == ", side);
		if(side == Location.left){
			writeln("population left wing part ctrl point geometry");
			foreach(wp_idx, wing_part; wing_geometry.wing_parts){
				foreach(crd_idx;0.._chordwise_nodes){
					foreach(sp_idx; 0.._spanwise_chunks){
					//writeln("going_left_wing_part");
						wing_part.ctrl_chunks[crd_idx*_spanwise_chunks + sp_idx].ctrl_pt_y[] = -span_ctrl_pt[sp_idx*chunk_size..sp_idx*chunk_size + chunk_size]*span[];
						wing_part.ctrl_chunks[crd_idx*_spanwise_chunks + sp_idx].ctrl_pt_x[] = chord_ctrl_pt[crd_idx]*wing_part.chunks[sp_idx].chord[] - wing_part.ctrl_chunks[crd_idx*_spanwise_chunks + sp_idx].ctrl_pt_y[]*tan(wing_part.le_sweep_angle*PI/180);
						wing_part.ctrl_chunks[crd_idx*_spanwise_chunks + sp_idx].camber[] = _camber;
						wing_part.ctrl_chunks[crd_idx*_spanwise_chunks + sp_idx].ctrl_pt_z[] = 0.0;
					}
				}
			}
		}
		else if(side == Location.right){
			writeln("population left wing part ctrl point geometry");
			foreach(wp_idx, wing_part; wing_geometry.wing_parts){
				foreach(crd_idx;0.._chordwise_nodes){
					foreach(sp_idx; 0.._spanwise_chunks){
						//writeln("going_left_wing_part");
						wing_part.ctrl_chunks[crd_idx*_spanwise_chunks + sp_idx].ctrl_pt_y[] = span_ctrl_pt[sp_idx*chunk_size..sp_idx*chunk_size + chunk_size]*span[];
						wing_part.ctrl_chunks[crd_idx*_spanwise_chunks + sp_idx].ctrl_pt_x[] = chord_ctrl_pt[crd_idx]*wing_part.chunks[sp_idx].chord[] + wing_part.ctrl_chunks[crd_idx*_spanwise_chunks + sp_idx].ctrl_pt_y[]*tan(wing_part.le_sweep_angle*PI/180);
						wing_part.ctrl_chunks[crd_idx*_spanwise_chunks + sp_idx].camber[] = _camber;
						wing_part.ctrl_chunks[crd_idx*_spanwise_chunks + sp_idx].ctrl_pt_z[] = 0.0;
					}
				}
			}
		}
	}else{
		foreach(wp_idx, wing_part; wing_geometry.wing_parts){
			foreach(crd_idx;0.._chordwise_nodes){
				foreach(sp_idx; 0.._spanwise_chunks){
					if(wp_idx%2==0){
					//writeln("going_left_wing_part");
						wing_part.ctrl_chunks[crd_idx*_spanwise_chunks + sp_idx].ctrl_pt_y[] = -span_ctrl_pt[sp_idx*chunk_size..sp_idx*chunk_size + chunk_size]*span[];
						wing_part.ctrl_chunks[crd_idx*_spanwise_chunks + sp_idx].ctrl_pt_x[] = chord_ctrl_pt[crd_idx]*wing_part.chunks[sp_idx].chord[] - wing_part.ctrl_chunks[crd_idx*_spanwise_chunks + sp_idx].ctrl_pt_y[]*tan(wing_part.le_sweep_angle*PI/180);
						wing_part.ctrl_chunks[crd_idx*_spanwise_chunks + sp_idx].camber[] = _camber;
					}else{
					//writeln("going_right_wing_part");
						wing_part.ctrl_chunks[crd_idx*_spanwise_chunks + sp_idx].ctrl_pt_y[] = span_ctrl_pt[sp_idx*chunk_size..sp_idx*chunk_size + chunk_size]*span[];
						wing_part.ctrl_chunks[crd_idx*_spanwise_chunks + sp_idx].ctrl_pt_x[] = chord_ctrl_pt[crd_idx]*wing_part.chunks[sp_idx].chord[] + wing_part.ctrl_chunks[crd_idx*_spanwise_chunks + sp_idx].ctrl_pt_y[]*tan(wing_part.le_sweep_angle*PI/180);
						wing_part.ctrl_chunks[crd_idx*_spanwise_chunks + sp_idx].camber[] = _camber;
					}				
				
					wing_part.ctrl_chunks[crd_idx*_spanwise_chunks + sp_idx].ctrl_pt_z[] = 0.0;
				}
			}
		}
	}    
}


/*double[] generate_chordwise_votex_nodes(size_t chord_n_sections) {
	import std.algorithm : map;
	import std.array : array;
	import std.math : cos, PI;
	import std.range : iota, retro;
	// Spanwise Votex nodes
	immutable num_points = chord_n_sections%chunk_size == 0 ? chord_n_sections : chord_n_sections + (chunk_size - chord_n_sections%chunk_size);
	return iota(1.0,num_points + 2.0).map!((k){
		immutable theta = (2*k-1)*PI/(2.0*(num_points.to!double + 1.0));
		auto chord_vortex_pt = 0.5*(1 - cos(theta)).to!double;
		return chord_vortex_pt;
	}).array;
}*/

double[] generate_chordwise_control_points(size_t chord_n_sections) {
	import std.algorithm : map;
	import std.array : array;
	import std.math : cos, PI;
	import std.range : iota, retro;
	// Spanwise Votex nodes
	immutable num_points = chord_n_sections%chunk_size == 0 ? chord_n_sections : chord_n_sections + (chunk_size - chord_n_sections%chunk_size);
	return iota(1.0,num_points + 2.0).map!((i){
		immutable theta = i*PI/(num_points.to!double + 1.0);
		auto chord_ctrl_pt = 0.5*(1 - cos(theta)).to!double;
		return chord_ctrl_pt;
	}).array;
}

extern(C++) struct Frame {
	Mat4 local_matrix;
	Mat4 global_matrix;
	
	NullableRef!Frame parent;

	this(Vec3 axis, double angle, Vec3 translation, Frame* _parent) {
		parent = nullableRef(_parent);
	}

	void rotate(Vec3 axis, double angle) {

	}

	void translate(Vec3 translation) {

	}
}

template is_aircraft(A) {
	enum bool is_aircraft = {
		static if(isPointer!(A)) {
			return isInstanceOf!(AircraftT, PointerTarget!A);
		} else {
			return isInstanceOf!(AircraftT, A);
		}
	}();
}

extern(C++) alias Aircraft = AircraftT!(ArrayContainer.none);

extern(C++) struct AircraftT(ArrayContainer AC) {

	mixin ArrayDeclMixin!(AC, RotorGeometryT!AC, "rotors");
	mixin ArrayDeclMixin!(AC, WingGeometryT!AC, "wings");
	//RotorGeometryT!AC[] rotors;

	this(size_t num_rotors , size_t num_wings) {		
		mixin(array_ctor_mixin!(AC, "RotorGeometryT!AC", "rotors", "num_rotors"));
		//rotors = new RotorGeometryT!AC[num_rotors];
		mixin(array_ctor_mixin!(AC, "WingGeometryT!AC", "wings", "num_wings"));
	}

	ref typeof(this) opAssign(typeof(this) ac) {
		import std.stdio : writeln;
		debug writeln("Aircraft opAssign");
		this.rotors = ac.rotors;
		this.wings = ac.wings;
		return this;
	}

	ref typeof(this) opAssign(ref typeof(this) ac) {
		this.rotors = ac.rotors;
		this.wings = ac.wings;
		return this;
	}

	ref typeof(this) opAssign(typeof(this)* ac) {
		this.rotors = ac.rotors;
		this.wings = ac.wings;
		return this;
	}

}

template is_rotor_geometry(A) {
	enum bool is_rotor_geometry = {
		static if(isPointer!(A)) {
			return isInstanceOf!(RotorGeometryT, PointerTarget!A);
		} else {
			return isInstanceOf!(RotorGeometryT, A);
		}
	}();
}

alias RotorGeometry = RotorGeometryT!(ArrayContainer.none);

extern (C++) struct RotorGeometryT(ArrayContainer AC) {

	mixin ArrayDeclMixin!(AC, BladeGeometryT!(AC), "blades");

	Vec3 origin;
	double radius;
	double solidity;

	this(size_t num_blades, Vec3 origin, double radius, double solidity) {
		mixin(array_ctor_mixin!(AC, "BladeGeometryT!(AC)", "blades", "num_blades"));

		this.origin = origin;
		this.radius = radius;
		this.solidity = solidity;
	}

	ref typeof(this) opAssign(typeof(this) rotor) {
		this.blades = rotor.blades;
		this.origin = rotor.origin;
		this.radius = rotor.radius;
		this.solidity = rotor.solidity;
		return this;
	}

	ref typeof(this) opAssign(ref typeof(this) rotor) {
		this.blades = rotor.blades;
		this.origin = rotor.origin;
		this.radius = rotor.radius;
		this.solidity = rotor.solidity;
		return this;
	}

	ref typeof(this) opAssign(typeof(this)* rotor) {
		this.blades = rotor.blades;
		this.origin = rotor.origin;
		this.radius = rotor.radius;
		this.solidity = rotor.solidity;
		return this;
	}

}

template is_blade_geometry_chunk(A) {
	enum bool is_blade_geometry_chunk = {
		static if(isPointer!(A)) {
			return is(PointerTarget!A == BladeGeometryChunk);
		} else {
			return is(A == BladeGeometryChunk);
		}
	}();
}

extern (C++) struct BladeGeometryChunk {
	/++
	 +   Twist distribution
	 +/
	Chunk twist; 
	/++
	 +   Chord distribution
	 +/
	Chunk chord;
	/++
	 +   Radial distribution
	 +/
	Chunk r;
	/++
	 +  Radial sectional airfoil lift curve slope
	 +/
	Chunk C_l_alpha;
	/++
	 +  Radial sectional 0 lift angle of attack
	 +/
	Chunk alpha_0;
	/++
	 +	Blade quarter chord sweep angle
	 +/
	 Chunk sweep;
}

template is_blade_geometry(A) {
	enum bool is_blade_geometry = {
		static if(isPointer!(A)) {
			return isInstanceOf!(BladeGeometryT, PointerTarget!A);
		} else {
			return isInstanceOf!(BladeGeometryT, A);
		}
	}();
}

alias BladeGeometry = BladeGeometryT!(ArrayContainer.none);
extern (C++) struct BladeGeometryT(ArrayContainer AC) {
	mixin ArrayDeclMixin!(AC, BladeGeometryChunk, "chunks");

	double azimuth_offset;
	double average_chord;

	BladeAirfoil airfoil;

	this(size_t num_elements, double azimuth_offset, double average_chord, BladeAirfoil airfoil) {
		immutable actual_num_elements = num_elements%chunk_size == 0 ? num_elements : num_elements + (chunk_size - num_elements%chunk_size);
		//enforce(num_elements % chunk_size == 0, "Number of spanwise elements must be a multiple of the chunk size ("~chunk_size.to!string~")");
		immutable num_chunks = actual_num_elements/chunk_size;
		mixin(array_ctor_mixin!(AC, "BladeGeometryChunk", "chunks", "num_chunks"));

		this.airfoil = airfoil;

		this.azimuth_offset = azimuth_offset;
		this.average_chord = average_chord;
	}

	ref typeof(this) opAssign(typeof(this) blade) {
		this.airfoil = blade.airfoil;
		this.chunks = blade.chunks;
		this.azimuth_offset = blade.azimuth_offset;
		this.average_chord = blade.average_chord;
		return this;
	}

	ref typeof(this) opAssign(ref typeof(this) blade) {
		this.airfoil = blade.airfoil;
		this.chunks = blade.chunks;
		this.azimuth_offset = blade.azimuth_offset;
		this.average_chord = blade.average_chord;
		return this;
	}

	ref typeof(this) opAssign(typeof(this)* blade) {
		this.airfoil = blade.airfoil;
		this.chunks = blade.chunks;
		this.azimuth_offset = blade.azimuth_offset;
		this.average_chord = blade.average_chord;
		return this;
	}
}


//Adding wing geometry

template is_wing_geometry(A) {
	enum bool is_wing_geometry = {
		static if(isPointer!(A)) {
			return isInstanceOf!(WingGeometryT, PointerTarget!A);
		} else {
			return isInstanceOf!(WingGeometryT, A);
		}
	}();
}

alias WingGeometry = WingGeometryT!(ArrayContainer.none);

extern (C++) struct WingGeometryT(ArrayContainer AC) {

	mixin ArrayDeclMixin!(AC, WingPartGeometryT!(AC), "wing_parts");

	Vec3 origin;
	double wing_span;

	this(size_t num_parts, Vec3 origin, double wing_span) {
		mixin(array_ctor_mixin!(AC, "WingPartGeometryT!(AC)", "wing_parts", "num_parts"));
		this.origin = origin;
	}

	ref typeof(this) opAssign(typeof(this) wing) {
		this.wing_parts = wing.wing_parts;
		this.origin = wing.origin;
		return this;
	}

	ref typeof(this) opAssign(ref typeof(this) wing) {
		this.wing_parts = wing.wing_parts;
		this.origin = wing.origin;
		return this;
	}

	ref typeof(this) opAssign(typeof(this)* wing) {
		this.wing_parts = wing.wing_parts;
		this.origin = wing.origin;
		return this;
	}

}

template is_wing_part_geometry_chunk(A) {
	enum bool is_wing_part_geometry_chunk = {
		static if(isPointer!(A)) {
			return is(PointerTarget!A == WingPartGeometryChunk);
		} else {
			return is(A == WingPartGeometryChunk);
		}
	}();
}

extern (C++) struct WingPartGeometryChunk {
	/++
	 +   Twist distribution
	 +/
	Chunk twist; 
	/++
	 +   Chord distribution
	 +/
	Chunk chord;
	//Chunk span_vortex_pt;

	//Chunk span_ctrl_pt;
	/++
	 +   Radial distribution
	 +/
	//Chunk r;
	/++
	 +  Radial sectional airfoil lift curve slope
	 +/
	//Chunk C_l_alpha;
	/++
	 +  Radial sectional 0 lift angle of attack
	 +/
	//Chunk alpha_0;
	/++
	 +	Blade quarter chord sweep angle
	 +/
	Chunk sweep;
}

extern (C++) struct WingPartCtrlPointChunk {
	
	Chunk ctrl_pt_x;

	Chunk ctrl_pt_y;

	Chunk ctrl_pt_z;

	Chunk camber;
}

template is_wing_part_geometry(A) {
	enum bool is_wing_geometry = {
		static if(isPointer!(A)) {
			return isInstanceOf!(WingPartGeometryT, PointerTarget!A);
		} else {
			return isInstanceOf!(WingPartGeometryT, A);
		}
	}();
}

enum Location : string{
		right = "right",
		left = "left",
		mid = "mid"
	}

alias WingPartGeometry = WingPartGeometryT!(ArrayContainer.none);
struct WingPartGeometryT(ArrayContainer AC) {
	mixin ArrayDeclMixin!(AC, WingPartGeometryChunk, "chunks");
	mixin ArrayDeclMixin!(AC, WingPartCtrlPointChunk, "ctrl_chunks");

	Vec3 wing_root_origin;
	double average_chord;
	double wing_root_chord;
	double wing_tip_chord;
	double le_sweep_angle;
	double te_sweep_angle;
	double wing_span;
	Location loc;

	//BladeAirfoil airfoil;

	this(size_t span_elements, size_t chordwise_nodes, Vec3 wing_root_origin, double average_chord, double wing_root_chord, double wing_tip_chord, double le_sweep_angle, double te_sweep_angle, double wing_span, Location loc) {
		immutable actual_num_span_elements = span_elements%chunk_size == 0 ? span_elements : span_elements + (chunk_size - span_elements%chunk_size);
		//enforce(num_elements % chunk_size == 0, "Number of spanwise elements must be a multiple of the chunk size ("~chunk_size.to!string~")");
		immutable num_span_chunks = actual_num_span_elements/chunk_size;
		//immutable actual_ctrl_pt = ((actual_num_span_elements)*(chord_elements))%chunk_size == 0 ? (span_elements-1)*(chord_elements) : (span_elements-1)*(chord_elements) + (chunk_size - (span_elements-1)*(chord_elements)%chunk_size);
		immutable ctrl_pt_chunks = num_span_chunks*chordwise_nodes;
		mixin(array_ctor_mixin!(AC, "WingPartGeometryChunk", "chunks", "num_span_chunks"));
		mixin(array_ctor_mixin!(AC, "WingPartCtrlPointChunk", "ctrl_chunks", "ctrl_pt_chunks"));

		//this.airfoil = airfoil;

		this.wing_root_origin = wing_root_origin;
		this.average_chord = average_chord;
		this.wing_root_chord = wing_root_chord;
		this.wing_tip_chord = wing_tip_chord;
		this.le_sweep_angle = le_sweep_angle;
		this.te_sweep_angle = te_sweep_angle;
		this.wing_span = wing_span;
		this.loc = loc;
	}

	ref typeof(this) opAssign(typeof(this) wing_part) {
		//this.airfoil = wing_part.airfoil;
		this.chunks = wing_part.chunks;
		this.ctrl_chunks = wing_part.ctrl_chunks;
		this.wing_root_origin = wing_part.wing_root_origin;
		this.average_chord = wing_part.average_chord;
		this.wing_root_chord = wing_part.wing_root_chord;
		this.wing_tip_chord = wing_part.wing_tip_chord;
		this.le_sweep_angle = wing_part.le_sweep_angle;
		this.te_sweep_angle = wing_part.te_sweep_angle;
		this.wing_span = wing_part.wing_span;
		this.loc = loc;
		return this;
	}

	ref typeof(this) opAssign(ref typeof(this) wing_part) {
		//this.airfoil = wing_part.airfoil;
		this.chunks = wing_part.chunks;
		this.ctrl_chunks = wing_part.ctrl_chunks;
		this.wing_root_origin = wing_part.wing_root_origin;
		this.average_chord = wing_part.average_chord;
		this.wing_root_chord = wing_part.wing_root_chord;
		this.wing_tip_chord = wing_part.wing_tip_chord;
		this.le_sweep_angle = wing_part.le_sweep_angle;
		this.te_sweep_angle = wing_part.te_sweep_angle;
		this.wing_span = wing_part.wing_span;
		this.loc = loc;
		return this;
	}

	ref typeof(this) opAssign(typeof(this)* wing_part) {
		//this.airfoil = wing_part.airfoil;
		this.chunks = wing_part.chunks;
		this.ctrl_chunks = wing_part.ctrl_chunks;
		this.wing_root_origin = wing_part.wing_root_origin;
		this.average_chord = wing_part.average_chord;
		this.wing_root_chord = wing_part.wing_root_chord;
		this.wing_tip_chord = wing_part.wing_tip_chord;
		this.le_sweep_angle = wing_part.le_sweep_angle;
		this.te_sweep_angle = wing_part.te_sweep_angle;
		this.wing_span = wing_part.wing_span;
		this.loc = loc;
		return this;
	}
}

void set_geometry_array(string value, ArrayContainer AC)(ref BladeGeometryT!AC blade, double[] data) {

	foreach(c_idx, ref chunk; blade.chunks) {
		immutable out_start_idx = c_idx*chunk_size;

		immutable remaining = data.length - out_start_idx;
		
		immutable out_end_idx = remaining > chunk_size ? (c_idx + 1)*chunk_size : out_start_idx + remaining;
		immutable in_end_idx = remaining > chunk_size ? chunk_size : remaining;

		mixin("chunk."~value~"[0..in_end_idx] = data[out_start_idx..out_end_idx];");
	}
}

void set_geometry_array(string value, ArrayContainer AC)(ref WingPartGeometryT!AC wing_part, double[] data) {

	foreach(c_idx, ref chunk; wing_part.chunks) {
		immutable out_start_idx = c_idx*chunk_size;

		immutable remaining = data.length - out_start_idx;
		
		immutable out_end_idx = remaining > chunk_size ? (c_idx + 1)*chunk_size : out_start_idx + remaining;
		immutable in_end_idx = remaining > chunk_size ? chunk_size : remaining;

		mixin("chunk."~value~"[0..in_end_idx] = data[out_start_idx..out_end_idx];");
	}
}