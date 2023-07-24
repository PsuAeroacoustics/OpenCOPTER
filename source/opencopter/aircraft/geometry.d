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
	//RotorGeometryT!AC[] rotors;

	this(size_t num_rotors) {
		mixin(array_ctor_mixin!(AC, "RotorGeometryT!AC", "rotors", "num_rotors"));
		//rotors = new RotorGeometryT!AC[num_rotors];
	}

	ref typeof(this) opAssign(typeof(this) ac) {
		import std.stdio : writeln;
		debug writeln("Aircraft opAssign");
		this.rotors = ac.rotors;
		return this;
	}

	ref typeof(this) opAssign(ref typeof(this) ac) {
		this.rotors = ac.rotors;
		return this;
	}

	ref typeof(this) opAssign(typeof(this)* ac) {
		this.rotors = ac.rotors;
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

void set_geometry_array(string value, ArrayContainer AC)(ref BladeGeometryT!AC blade, double[] data) {

	foreach(c_idx, ref chunk; blade.chunks) {
		immutable out_start_idx = c_idx*chunk_size;

		immutable remaining = data.length - out_start_idx;
		
		immutable out_end_idx = remaining > chunk_size ? (c_idx + 1)*chunk_size : out_start_idx + remaining;
		immutable in_end_idx = remaining > chunk_size ? chunk_size : remaining;

		mixin("chunk."~value~"[0..in_end_idx] = data[out_start_idx..out_end_idx];");
	}
}
