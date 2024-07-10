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
import std.stdio : writeln;

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

		translate(translation);
		rotate(axis, angle);
	}

	Vec3 global_position() {
		nop;
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
		M[1, 0] = xy*one_cs+z*sn; M[1, 1] = cs+yy*one_cs;   M[1, 2]  = yz-x*sn;		 
		M[2, 0] = xz*one_cs-y*sn; M[2, 1] = yz*one_cs+x*sn; M[2, 2]  = cs+zz*one_cs;	 
		
		temp_mat = temp_mat * local_matrix.extract_rotation_matrix();
		
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
		M[1, 0] = xy*one_cs+z*sn; M[1, 1] = cs+yy*one_cs;   M[1, 2]  = yz-x*sn;		 
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

	void translate(Vec3 translation) {
		nop;
		local_matrix[0, 3] += translation[0];
		local_matrix[1, 3] += translation[1];
		local_matrix[2, 3] += translation[2];
	}

	void update(ref Mat4 parent_global_mat) {
		global_matrix = parent_global_mat*local_matrix;
		inverse_global_matrix = global_matrix.inverse.get;

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

	Frame* root_frame;

	this(size_t num_rotors) {
		mixin(array_ctor_mixin!(AC, "RotorGeometryT!AC", "rotors", "num_rotors"));

		root_frame = new Frame(Vec3(0, 0, 1), PI, Vec3(0, 0, 0), null, "aircraft", "connection");
	}

	ref typeof(this) opAssign(typeof(this) ac) {
		import std.stdio : writeln;
		debug writeln("Aircraft opAssign");
		this.rotors = ac.rotors;
		this.root_frame = ac.root_frame;
		return this;
	}

	ref typeof(this) opAssign(ref typeof(this) ac) {
		this.rotors = ac.rotors;
		this.root_frame = ac.root_frame;
		return this;
	}

	ref typeof(this) opAssign(typeof(this)* ac) {
		this.rotors = ac.rotors;
		this.root_frame = ac.root_frame;
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

void compute_blade_vectors(BG)(ref BG blade) {

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

alias compute_blade_vectors_test = compute_blade_vectors!(BladeGeometry);

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
