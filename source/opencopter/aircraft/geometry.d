module opencopter.aircraft.geometry;

import opencopter.aircraft;
import opencopter.config;
import opencopter.math;
import opencopter.memory;
import opencopter.airfoilmodels;

import numd.linearalgebra.matrix;

import std.conv : to;
import std.exception : enforce;
import std.math;
import std.range;
import std.traits;
import std.typecons;
import std.stdio;


double[] generate_radius_points(size_t n_sections, double root_cutout = 0.0) {
	import std.algorithm : map;
	import std.array : array;
	import std.math : cos, PI;
	import std.range : iota, retro;

	immutable num_points = n_sections%chunk_size == 0 ? n_sections : n_sections + (chunk_size - n_sections%chunk_size);
    return iota(1.0, num_points + 1.0).map!((n) {
    	immutable psi = n*PI/(num_points.to!double + 1.0);
    	auto r = (1.0 - root_cutout)*0.5*(cos(psi) + 1.0).to!double + root_cutout;
    	return r;
    }).retro.array;
}

// double[] generate_radius_points_half_cos(size_t n_sections, double root_cutout = 0.0) {
// 	import std.algorithm : map;
// 	import std.array : array;
// 	import std.math : cos, PI;
// 	import std.range : iota, retro;

// 	immutable num_points = n_sections%chunk_size == 0 ? n_sections : n_sections + (chunk_size - n_sections%chunk_size);
//     //return iota(1.0, num_points + 1.0).map!((n) {
// 	return iota(0, 0.5*PI, 0.5*PI/num_points.to!double).map!((psi) {
//     	//immutable psi = n*PI/(num_points.to!double + 1.0);
//     	auto r = (1.0 - root_cutout)*(cos(psi) + 1.0).to!double + root_cutout;
//     	return r;
//     }).retro.array;
// }

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
	size_t acutal_span_nodes = _spanwise_nodes%chunk_size == 0 ? _spanwise_nodes : _spanwise_nodes + (chunk_size - _spanwise_nodes%chunk_size);
	size_t _spanwise_chunks = acutal_span_nodes/chunk_size;
    auto span_ctrl_pt = generate_spanwise_control_points(acutal_span_nodes);
	auto chord_ctrl_pt = generate_chordwise_control_points(_chordwise_nodes);
    //auto span_ctrl_pt_left = generate_spanwise_control_points(_spanwise_chunks*chunk_size);
    //span_vr_nodes_left = span_vr_nodes_left.retro;
    //auto chord_vr_nodes = generate_chordwise_votex_nodes(_chordwise_nodes);

	foreach(wp_idx, wing_part; wing_geometry.wing_parts){
		// A lone wing part lies on the side its location names. With several
		// parts the even ones are on the left (negative y) and the odd ones on
		// the right, as set_wing_vortex_geometry places them.
		double side;
		if(wing_geometry.wing_parts.length == 1){
			string loc = wing_part.loc;
			writeln("location == ", loc);
			if(loc == Location.left){
				writeln("population left wing part ctrl point geometry");
				side = -1.0;
			}else if(loc == Location.right){
				writeln("populating right wing part ctrl point geometry");
				side = 1.0;
			}else{
				continue;
			}
		}else{
			side = wp_idx%2 == 0 ? -1.0 : 1.0;
		}

		// Each part uses its own half span. Positive LE sweep is aft on both
		// sides.
		immutable double span = 2*wing_part.wing_span;
		immutable tan_le_sweep = tan(wing_part.le_sweep_angle*PI/180);

		foreach(crd_idx;0.._chordwise_nodes){
			foreach(sp_idx; 0.._spanwise_chunks){
				immutable Chunk abs_y = span_ctrl_pt[sp_idx*chunk_size..sp_idx*chunk_size + chunk_size]*span;
				wing_part.ctrl_chunks[crd_idx*_spanwise_chunks + sp_idx].ctrl_pt_y[] = side*abs_y[];
				wing_part.ctrl_chunks[crd_idx*_spanwise_chunks + sp_idx].ctrl_pt_x[] = chord_ctrl_pt[crd_idx]*wing_part.chunks[sp_idx].chord[] + abs_y[]*tan_le_sweep;
				wing_part.ctrl_chunks[crd_idx*_spanwise_chunks + sp_idx].camber[] = _camber;
				wing_part.ctrl_chunks[crd_idx*_spanwise_chunks + sp_idx].ctrl_pt_z[] = 0.0;
			}
		}
	}
}

double[] generate_chordwise_control_points(size_t chord_n_sections) {
	import std.algorithm : map;
	import std.array : array;
	import std.math : cos, PI;
	import std.range : iota;
	// Chordwise control points of Lan's quasi vortex lattice at
	// x/c = (1 - cos(theta_i))/2, theta_i = i*PI/N for i = 1 .. N, the
	// stations VortexLatticeT's influence matrix is built for. The last one
	// is on the trailing edge. Chordwise control points are not stored in
	// chunks, so N is not rounded up to a multiple of chunk_size.
	return iota(1, chord_n_sections + 1).map!((i){
		immutable theta = i*PI/(chord_n_sections.to!double);
		auto chord_ctrl_pt = 0.5*(1 - cos(theta)).to!double;
		return chord_ctrl_pt;
	}).array;
}

Mat3 extract_rotation_matrix(ref Mat4 mat) {
	Mat3 rot_mat;

	rot_mat[0, 0] = mat[0, 0];
	rot_mat[0, 1] = mat[0, 1];
	rot_mat[0, 2] = mat[0, 2];

	rot_mat[1, 0] = mat[1, 0];
	rot_mat[1, 1] = mat[1, 1];
	rot_mat[1, 2] = mat[1, 2];

	rot_mat[2, 0] = mat[2, 0];
	rot_mat[2, 1] = mat[2, 1];
	rot_mat[2, 2] = mat[2, 2];

	return rot_mat;
}

void set_rotation_matrix(ref Mat4 mat4, ref Mat3 rot_mat) {
	mat4[0, 0] = rot_mat[0, 0];
	mat4[0, 1] = rot_mat[0, 1];
	mat4[0, 2] = rot_mat[0, 2];

	mat4[1, 0] = rot_mat[1, 0];
	mat4[1, 1] = rot_mat[1, 1];
	mat4[1, 2] = rot_mat[1, 2];

	mat4[2, 0] = rot_mat[2, 0];
	mat4[2, 1] = rot_mat[2, 1];
	mat4[2, 2] = rot_mat[2, 2];
}

Mat3 build_rotation_matrix(Vec3 axis, double angle) {
	Mat3 temp_mat;
	double cs = cos(angle);
	double sn = sin(angle);

	double x = axis[0];
	double y = axis[1];
	double z = axis[2];

	double mag = sqrt((x*x)+(y*y)+(z*z));
	x = x/mag;
	y = y/mag;
	z = z/mag;
	
	double xx = x*x;
	double xy = x*y;
	double xz = x*z;
	double yy = y*y;
	double yz = y*z;
	double zz = z*z;
	double one_cs = 1 - cs;
	
	alias M = temp_mat;
	M[0, 0] = cs + xx*one_cs; M[0, 1] = xy*one_cs-z*sn; M[0, 2]  = xz*one_cs+y*sn;
	M[1, 0] = xy*one_cs+z*sn; M[1, 1] = cs+yy*one_cs;   M[1, 2]  = yz*one_cs-x*sn;		 
	M[2, 0] = xz*one_cs-y*sn; M[2, 1] = yz*one_cs+x*sn; M[2, 2]  = cs+zz*one_cs;	 

	return temp_mat;
}

void print_frame(F)(F frame, int depth = 0) {
	import std.stdio : writeln;
	import std.range : repeat;

	writeln("\t".repeat(depth).join, " ", frame.name, ": ", frame.global_matrix);//, "\t", frame.local_matrix);

	foreach(ref child; frame.children) {
		print_frame(child, depth + 1);
	}
}

enum FrameType : int {
	aircraft,
	connection,
	rotor,
	blade,
	wing
}

struct Frame {

	Mat4 local_matrix;
	Mat4 global_matrix;
	Mat4 inverse_global_matrix;
	string name;
	FrameType frame_type;

	Frame* parent;
	Frame*[] children;

	Vec3 axis;
	double angle;

	this(Vec3 _axis, double _angle, Vec3 translation, Frame* _parent, string _name, string _frame_type) {
		parent = _parent;
		name = _name;
		axis = _axis;
		angle = _angle;

		frame_type = _frame_type.to!FrameType;
		
		local_matrix = Mat4.identity();
		global_matrix = Mat4.identity();
		inverse_global_matrix = Mat4.identity();

		translate(translation);
		rotate(axis, angle);
	}

	this(Frame* _frame, string _name) {
		this.local_matrix = _frame.local_matrix;
		this.global_matrix = _frame.global_matrix;
		this.name = _name;
		this.frame_type = _frame.frame_type;
		this.parent = _frame.parent;
		this.children = _frame.children.dup;
		this.axis = _frame.axis;
		this.angle = _frame.angle;
	}

	this(Vec3 _axis, double _angle, Vec3 translation, Frame* _parent, string _name, FrameType _frame_type) {
		parent = _parent;
		name = _name;
		axis = _axis;
		angle = _angle;

		frame_type = _frame_type;
		
		local_matrix = Mat4.identity();
		global_matrix = Mat4.identity();
		inverse_global_matrix = Mat4.identity();

		translate(translation);
		rotate(axis, angle);
	}

	Vec3 global_position() {
		nop;
		// The translation of the global matrix is this frame's origin in global
		// coordinates (the inverse holds the global origin in local coordinates).
		auto pos = Vec3(global_matrix[0, 3], global_matrix[1, 3], global_matrix[2, 3]);
		return pos;
	}

	Vec3 local_position() {
		nop;
		auto pos = Vec3(local_matrix[0, 3], local_matrix[1, 3], local_matrix[2, 3]);
		return pos;
	}

	string get_frame_type() {
		nop;
		return frame_type.to!string;
	}

	void set_frame_type(string _frame_type) {
		nop;
		frame_type = _frame_type.to!FrameType;
	}

	void rotate(Vec3 axis, double angle) {
		Mat3 temp_mat;
		double cs = cos(angle);
		double sn = sin(angle);

		double x = axis[0];
		double y = axis[1];
		double z = axis[2];

		double mag = sqrt((x*x)+(y*y)+(z*z));
		x = x/mag;
		y = y/mag;
		z = z/mag;
		
		double xx = x*x;
		double xy = x*y;
		double xz = x*z;
		double yy = y*y;
		double yz = y*z;
		double zz = z*z;
		double one_cs = 1 - cs;
		
		alias M = temp_mat;
		M[0, 0] = cs + xx*one_cs; M[0, 1] = xy*one_cs-z*sn; M[0, 2]  = xz*one_cs+y*sn;
		M[1, 0] = xy*one_cs+z*sn; M[1, 1] = cs+yy*one_cs;   M[1, 2]  = yz*one_cs-x*sn;		 
		M[2, 0] = xz*one_cs-y*sn; M[2, 1] = yz*one_cs+x*sn; M[2, 2]  = cs+zz*one_cs;	 
		
		temp_mat = local_matrix.extract_rotation_matrix() * temp_mat;
		
		local_matrix.set_rotation_matrix(temp_mat);
	}

	void set_rotation(Vec3 _axis, double _angle) {

		axis = _axis;
		angle = _angle;

		Mat3 temp_mat;
		double cs = cos(angle);
		double sn = sin(angle);

		double x = axis[0];
		double y = axis[1];
		double z = axis[2];

		double mag = sqrt((x*x)+(y*y)+(z*z));
		x = x/mag;
		y = y/mag;
		z = z/mag;
		
		double xx = x*x;
		double xy = x*y;
		double xz = x*z;
		double yy = y*y;
		double yz = y*z;
		double zz = z*z;
		double one_cs = 1 - cs;
		
		alias M = temp_mat;
		M[0, 0] = cs + xx*one_cs; M[0, 1] = xy*one_cs-z*sn; M[0, 2]  = xz*one_cs+y*sn;
		M[1, 0] = xy*one_cs+z*sn; M[1, 1] = cs+yy*one_cs;   M[1, 2]  = yz*one_cs-x*sn;		 
		M[2, 0] = xz*one_cs-y*sn; M[2, 1] = yz*one_cs+x*sn; M[2, 2]  = cs+zz*one_cs;	 
		
		local_matrix.set_rotation_matrix(temp_mat);
	}

	Vec3 local_rotation_axis() {
		nop;
		return axis;
	}

	double local_rotation_angle() {
		nop;
		return angle;
	}

	Mat4 get_local_mat(){
		nop;
		return local_matrix;
	}

	Mat4 get_global_mat(){
		nop;
		return global_matrix;
	}

	Vec4 get_origin_in_frame_coordinate(Vec4 origin){
		auto local_origin = -local_matrix*origin;
		return local_origin;
	}

	void translate(Vec3 translation) {
		nop;
		local_matrix[0, 3] += translation[0];
		local_matrix[1, 3] += translation[1];
		local_matrix[2, 3] += translation[2];
	}

	void update(ref Mat4 parent_global_mat) {
		global_matrix = parent_global_mat * local_matrix;
		auto rotation_matrix = extract_rotation_matrix(global_matrix);
		auto inverse_rotation_matrix = rotation_matrix.transpose;
		set_rotation_matrix(inverse_global_matrix, inverse_rotation_matrix);
		auto d_vec = Vec3(global_matrix[0, 3], global_matrix[1, 3], global_matrix[2, 3]);
		auto invese_d_vec = -inverse_rotation_matrix*d_vec;
		//writeln("inverse d vec: ", "\tx = ", invese_d_vec[0], "\ty = ", invese_d_vec[1], "\tz = ", invese_d_vec[2]);
		inverse_global_matrix[0, 3] = invese_d_vec[0];
		inverse_global_matrix[1, 3] = invese_d_vec[1];
		inverse_global_matrix[2, 3] = invese_d_vec[2];
		inverse_global_matrix[3, 3] = 1.0;
		
		foreach(ref child; children){
			child.update(global_matrix);
		}
	}

	void update(Mat4 parent_global_mat) {
		update(parent_global_mat);
	}

	void update(Mat4* parent_global_mat) {
		update(*parent_global_mat);
	}
	
	auto global_to_local(Vec3 global_pos) {
		nop;
		auto res = inverse_global_matrix*Vec4(global_pos[0], global_pos[1], global_pos[2], 1.0);
		return Vec3(res[0], res[1], res[2]);
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

	Frame* root_frame;

	this(size_t num_rotors , size_t num_wings) {		
		mixin(array_ctor_mixin!(AC, "RotorGeometryT!AC", "rotors", "num_rotors"));

		root_frame = new Frame(Vec3(0, 0, 1), 0.0, Vec3(0, 0, 0), null, "aircraft", "connection");
		mixin(array_ctor_mixin!(AC, "WingGeometryT!AC", "wings", "num_wings"));
	}

	ref typeof(this) opAssign(typeof(this) ac) {
		import std.stdio : writeln;
		debug writeln("Aircraft opAssign");
		this.rotors = ac.rotors;
		this.root_frame = ac.root_frame;
		this.wings = ac.wings;
		return this;
	}

	ref typeof(this) opAssign(ref typeof(this) ac) {
		this.rotors = ac.rotors;
		this.root_frame = ac.root_frame;
		this.wings = ac.wings;
		return this;
	}

	ref typeof(this) opAssign(typeof(this)* ac) {
		this.rotors = ac.rotors;
		this.root_frame = ac.root_frame;
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

	Frame* frame;

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
		this.frame = rotor.frame;
		return this;
	}

	ref typeof(this) opAssign(ref typeof(this) rotor) {
		this.blades = rotor.blades;
		this.origin = rotor.origin;
		this.radius = rotor.radius;
		this.solidity = rotor.solidity;
		this.frame = rotor.frame;
		return this;
	}

	ref typeof(this) opAssign(typeof(this)* rotor) {
		this.blades = rotor.blades;
		this.origin = rotor.origin;
		this.radius = rotor.radius;
		this.solidity = rotor.solidity;
		this.frame = rotor.frame;
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

	/++
	 +	Blade local normalized x offset;
	 +/
	Chunk xi;

	/++
	 +	Blade local normalized x offset derivative;
	 +/
	Chunk xi_p;

	Chunk thickness;
	
	Vector!(4, Chunk) af_norm;
}

Matrix!(r, c, double) extract_single_mat(size_t r, size_t c)(auto ref Matrix!(r, c, Chunk) mat, size_t chunk_idx) {
	Matrix!(r, c, double) ret_mat;

	for(size_t i = 0; i < r; i++) {
		for(size_t j = 0; j < c; j++) {
			ret_mat[i, j] = mat[i, j][chunk_idx];
		}
	}

	return ret_mat;
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

	/++
	 +	True length of the blade accounting for root cutout.
	 +/
	double blade_length;

	BladeAirfoil airfoil;

	Frame* frame;

	double r_c;

	this(size_t num_elements, double azimuth_offset, double average_chord, BladeAirfoil airfoil, double r_c) {
		immutable actual_num_elements = num_elements%chunk_size == 0 ? num_elements : num_elements + (chunk_size - num_elements%chunk_size);

		immutable num_chunks = actual_num_elements/chunk_size;
		mixin(array_ctor_mixin!(AC, "BladeGeometryChunk", "chunks", "num_chunks"));

		this.airfoil = airfoil;

		this.azimuth_offset = azimuth_offset;
		this.average_chord = average_chord;
		this.r_c = r_c;
	}

	this(ref typeof(this) blade) {
		this.airfoil = blade.airfoil;
		this.chunks = blade.chunks;
		this.azimuth_offset = blade.azimuth_offset;
		this.average_chord = blade.average_chord;
		this.blade_length = blade.blade_length;
		this.frame = blade.frame;
		this.r_c = blade.r_c;
	}

	this(typeof(this)* blade) {
		this.airfoil = blade.airfoil;
		this.chunks = blade.chunks;
		this.azimuth_offset = blade.azimuth_offset;
		this.average_chord = blade.average_chord;
		this.blade_length = blade.blade_length;
		this.frame = blade.frame;
		this.r_c = blade.r_c;
	}

	ref typeof(this) opAssign(typeof(this) blade) {
		this.airfoil = blade.airfoil;
		this.chunks = blade.chunks;
		this.azimuth_offset = blade.azimuth_offset;
		this.average_chord = blade.average_chord;
		this.blade_length = blade.blade_length;
		this.frame = blade.frame;
		this.r_c = blade.r_c;
		return this;
	}

	ref typeof(this) opAssign(ref typeof(this) blade) {
		this.airfoil = blade.airfoil;
		this.chunks = blade.chunks;
		this.azimuth_offset = blade.azimuth_offset;
		this.average_chord = blade.average_chord;
		this.blade_length = blade.blade_length;
		this.frame = blade.frame;
		this.r_c = blade.r_c;
		return this;
	}

	ref typeof(this) opAssign(typeof(this)* blade) {
		this.airfoil = blade.airfoil;
		this.chunks = blade.chunks;
		this.azimuth_offset = blade.azimuth_offset;
		this.average_chord = blade.average_chord;
		this.blade_length = blade.blade_length;
		this.frame = blade.frame;
		this.r_c = blade.r_c;
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

	Frame* frame;

	this(size_t num_parts, Vec3 origin, double wing_span) {
		mixin(array_ctor_mixin!(AC, "WingPartGeometryT!(AC)", "wing_parts", "num_parts"));
		this.origin = origin;
		this.wing_span = wing_span;
	}

	ref typeof(this) opAssign(typeof(this) wing) {
		this.wing_parts = wing.wing_parts;
		this.origin = wing.origin;
		this.wing_span = wing.wing_span;
		this.frame = wing.frame;
		return this;
	}

	ref typeof(this) opAssign(ref typeof(this) wing) {
		this.wing_parts = wing.wing_parts;
		this.origin = wing.origin;
		this.wing_span = wing.wing_span;
		this.frame = wing.frame;
		return this;
	}

	ref typeof(this) opAssign(typeof(this)* wing) {
		this.wing_parts = wing.wing_parts;
		this.origin = wing.origin;
		this.wing_span = wing.wing_span;
		this.frame = wing.frame;
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
	/++
	 +   span_locations
	 +/
	Chunk y_span;
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

	Frame* frame;
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
		this.loc = wing_part.loc;
		//this.frame = wing_part.frame;
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
		this.loc = wing_part.loc;
		//this.frame = wing_part.frame;
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
		this.loc = wing_part.loc;
		//this.frame = wing_part.frame;
		return this;
	}
}

void compute_blade_vectors(BG)(auto ref BG blade) {

	auto twist_array = blade.get_geometry_array!"twist";
	auto sweep_array = blade.get_geometry_array!"sweep";
	auto r_array = blade.get_geometry_array!"r";
	auto xi_array = blade.get_geometry_array!"xi";

	Frame af_frame;

	Vec4[] af_norms = new Vec4[twist_array.length];

	auto af_vec = Vec4(0);
	af_vec[1] = 1;

	auto twist_axis = Vec3(1, 0, 0);
	auto sweep_axis = Vec3(0, 0, 1);

	foreach(ref element; zip(sweep_array, twist_array, r_array, xi_array).enumerate()) {
		auto idx = element[0];
		auto sweep = element[1][0];
		auto twist = element[1][1];
		auto r = element[1][2];
		auto xi = element[1][3];

		af_frame.local_matrix = Mat4.identity();
		af_frame.global_matrix = Mat4.identity();
		af_frame.children ~= new Frame();
		af_frame.children[0].local_matrix = Mat4.identity();
		af_frame.children[0].global_matrix = Mat4.identity();

		immutable af_pos = Vec3(r, xi, 0);
		af_frame.translate(af_pos);
		af_frame.rotate(sweep_axis, sweep);
		af_frame.children[0].rotate(twist_axis, twist);
		af_frame.update(Mat4.identity);

		af_frame.rotate(sweep_axis, -0.5*PI);
		af_frame.update(Mat4.identity);

		af_norms[idx] = (af_frame.global_matrix*af_vec).normalize();
	}

	blade.set_geometry_array!"af_norm"(af_norms);
}


double[] sweep_from_quarter_chord(double[] r, double[] xi) {
	double[] sweep = new double[xi.length];

	double rise = 0;
	double run = 0;
	foreach(idx; 0..xi.length) {
		if(idx == 0) {
			rise = xi[idx + 1] - xi[idx];
			run = r[idx + 1] - r[idx];
		} else if(idx == xi.length - 1) {
			rise = xi[idx] - xi[idx - 1];
			run = r[idx] - r[idx - 1];
		} else {
			rise = xi[idx + 1] - xi[idx - 1];
			run = r[idx + 1] - r[idx - 1];
		}

		static import std.math;
		sweep[idx] = std.math.atan(-rise/run);
	}

	return sweep;
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

void set_geometry_array(string value, ArrayContainer AC)(ref BladeGeometryT!AC blade, Vec4[] data) {

	foreach(c_idx, ref chunk; blade.chunks) {
		immutable out_start_idx = c_idx*chunk_size;

		immutable remaining = data.length - out_start_idx;
		
		immutable out_end_idx = remaining > chunk_size ? (c_idx + 1)*chunk_size : out_start_idx + remaining;

		mixin("chunk."~value~" = data[out_start_idx..out_end_idx];");
	}
}

void set_geometry_array(string value, ArrayContainer AC)(ref BladeGeometryT!AC blade, Mat4[] data) {

	foreach(c_idx, ref chunk; blade.chunks) {
		immutable out_start_idx = c_idx*chunk_size;

		immutable remaining = data.length - out_start_idx;
		
		immutable out_end_idx = remaining > chunk_size ? (c_idx + 1)*chunk_size : out_start_idx + remaining;

		mixin("chunk."~value~" = data[out_start_idx..out_end_idx];");
	}
}

void set_geometry_array(string value, ArrayContainer AC)(BladeGeometryT!AC* blade, double[] data) {

	foreach(c_idx, ref chunk; blade.chunks) {
		immutable out_start_idx = c_idx*chunk_size;

		immutable remaining = data.length - out_start_idx;
		
		immutable out_end_idx = remaining > chunk_size ? (c_idx + 1)*chunk_size : out_start_idx + remaining;
		immutable in_end_idx = remaining > chunk_size ? chunk_size : remaining;

		mixin("chunk."~value~"[0..in_end_idx] = data[out_start_idx..out_end_idx];");
	}
}

void set_geometry_array(string value, ArrayContainer AC)(BladeGeometryT!AC* blade, Vec4[] data) {

	foreach(c_idx, ref chunk; blade.chunks) {
		immutable out_start_idx = c_idx*chunk_size;

		immutable remaining = data.length - out_start_idx;
		
		immutable out_end_idx = remaining > chunk_size ? (c_idx + 1)*chunk_size : out_start_idx + remaining;

		mixin("chunk."~value~" = data[out_start_idx..out_end_idx];");
	}
}

void set_geometry_array(string value, ArrayContainer AC)(BladeGeometryT!AC* blade, Mat4[] data) {

	foreach(c_idx, ref chunk; blade.chunks) {
		immutable out_start_idx = c_idx*chunk_size;

		immutable remaining = data.length - out_start_idx;
		
		immutable out_end_idx = remaining > chunk_size ? (c_idx + 1)*chunk_size : out_start_idx + remaining;

		mixin("chunk."~value~" = data[out_start_idx..out_end_idx];");
	}
}

double[] get_geometry_array(string value, ArrayContainer AC)(ref BladeGeometryT!AC blade) {
	immutable elements = blade.chunks.length*chunk_size;
	double[] geom_array = new double[elements];
	foreach(c_idx, ref chunk; blade.chunks) {
		immutable out_start_idx = c_idx*chunk_size;

		immutable remaining = elements - out_start_idx;
		
		immutable out_end_idx = remaining > chunk_size ? (c_idx + 1)*chunk_size : out_start_idx + remaining;
		immutable in_end_idx = remaining > chunk_size ? chunk_size : remaining;

		mixin("geom_array[out_start_idx..out_end_idx] = chunk."~value~"[0..in_end_idx];");
	}
	return geom_array;
}

double[] get_geometry_array(string value, ArrayContainer AC)(BladeGeometryT!AC* blade) {
	immutable elements = blade.chunks.length*chunk_size;
	double[] geom_array = new double[elements];
	foreach(c_idx, ref chunk; blade.chunks) {
		immutable out_start_idx = c_idx*chunk_size;

		immutable remaining = elements - out_start_idx;
		
		immutable out_end_idx = remaining > chunk_size ? (c_idx + 1)*chunk_size : out_start_idx + remaining;
		immutable in_end_idx = remaining > chunk_size ? chunk_size : remaining;

		mixin("geom_array[out_start_idx..out_end_idx] = chunk."~value~"[0..in_end_idx];");
	}
	return geom_array;
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


void set_geometry_array(string value, ArrayContainer AC)(WingPartGeometryT!AC* wing_part, double[] data) {

	foreach(c_idx, ref chunk; wing_part.chunks) {
		immutable out_start_idx = c_idx*chunk_size;

		immutable remaining = data.length - out_start_idx;
		
		immutable out_end_idx = remaining > chunk_size ? (c_idx + 1)*chunk_size : out_start_idx + remaining;
		immutable in_end_idx = remaining > chunk_size ? chunk_size : remaining;

		mixin("chunk."~value~"[0..in_end_idx] = data[out_start_idx..out_end_idx];");
	}
}

unittest {
	// Rotations about an oblique axis stay orthogonal (901722d).
	import std.math : abs, fmax;

	double orthogonality_error(Mat3 m) {
		double error = 0;
		foreach(i; 0..3) {
			foreach(j; 0..3) {
				double dot = 0;
				foreach(k; 0..3) {
					dot += m[k, i]*m[k, j];
				}
				error = fmax(error, abs(dot - (i == j ? 1.0 : 0.0)));
			}
		}
		return error;
	}

	foreach(axis; [Vec3(0, 0, 1), Vec3(0, 1, 1), Vec3(1, 2, 3)]) {
		assert(orthogonality_error(build_rotation_matrix(axis, 0.7)) < 1.0e-14);

		auto set = Frame(Vec3(1, 0, 0), 0.0, Vec3(0, 0, 0), null, "set", FrameType.connection);
		set.set_rotation(axis, 0.7);
		assert(orthogonality_error(extract_rotation_matrix(set.local_matrix)) < 1.0e-14);

		auto rotated = Frame(Vec3(1, 0, 0), 0.0, Vec3(0, 0, 0), null, "rotated", FrameType.connection);
		rotated.rotate(axis, 0.7);
		assert(orthogonality_error(extract_rotation_matrix(rotated.local_matrix)) < 1.0e-14);
	}
}

unittest {
	// global_position is the frame's origin in global coordinates (7d90fa4).
	import std.math : abs, PI;

	auto root = new Frame(Vec3(0, 0, 1), 0.0, Vec3(0, 0, 0), null, "root", FrameType.aircraft);
	auto hub = new Frame(Vec3(0, 0, 1), PI/2, Vec3(1, 2, 3), root, "hub", FrameType.connection);
	auto rotor = new Frame(Vec3(1, 0, 0), 0.0, Vec3(0.5, 0, 0), hub, "rotor", FrameType.rotor);
	root.children ~= hub;
	hub.children ~= rotor;
	root.update(Mat4.identity);
	auto position = rotor.global_position();
	auto origin = rotor.global_matrix*Vec4(0, 0, 0, 1);
	foreach(i; 0..3) {
		assert(abs(position[i] - origin[i]) < 1.0e-12);
	}
}

unittest {
	// Chordwise control points sit on Lan's stations, x/c = (1 - cos(i pi/N))/2,
	// with no rounding of N to whole chunks (7c8c820).
	import std.math : abs, cos, PI;

	foreach(N; [1UL, 2, 3, 4, 6, 8, 12]) {
		auto x = generate_chordwise_control_points(N);
		assert(x.length == N);
		foreach(i; 0..N) {
			assert(abs(x[i] - 0.5*(1.0 - cos((i + 1.0)*PI/N))) < 1.0e-15);
		}
	}
}

unittest {
	// Assigning a wing part keeps its side, and a wing keeps its span (5fd32ba).
	WingPartGeometry part;
	part = WingPartGeometry(8, 2, Vec3(0, 0, 0), 0.2, 0.2, 0.2, 0, 0, 1.0, Location.left);
	assert(part.loc == Location.left);
	assert(WingGeometry(1, Vec3(0, 0, 0), 3.5).wing_span == 3.5);
}
