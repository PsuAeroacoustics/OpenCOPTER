module opencopter.math.vectorarray;

import opencopter.math;

import numd.linearalgebra.matrix;

import std.algorithm;
import std.array;
import std.traits;

/++
 +/
struct VecChunkT(size_t dim, T) if(isStaticArray!T) {
	alias Vec = Vector!(dim, T);
	alias Mat = Matrix!(dim, dim, T);
	T[dim] data;
	alias data this;

	alias ThisType = typeof(this);

	/++
	+/
	this(Vector!(dim, T) init_data) {
		data[0][] = init_data[0];
		data[1][] = init_data[1];
		data[2][] = init_data[2];
	}

	this(immutable ThisType init_data) {
		data[0][] = init_data.data[0][];
		data[1][] = init_data.data[1][];
		data[2][] = init_data.data[2][];
	}

	@trusted @nogc auto get_sum() immutable {
		T[dim] d = data[].map!(a => a[].sum).staticArray!(T[dim]);
		return Vec(d);
	}

	@trusted @nogc ref auto get_vector(size_t i) immutable pure nothrow {
		T[dim] d = data[].map!(a => a[i]).staticArray!(T[dim]);
		return Vec(d);
	}

	@trusted @nogc auto opUnary(string op)() immutable pure nothrow
		if(op == "-")
	{
		auto ret = ThisType();

		static foreach(i; 0..dim) {
			ret[i][] = -data[i][];
		}

		return ret;
	}

	@trusted @nogc auto opBinary(string op)(ref ThisType rhs) immutable pure nothrow {
		auto ret = ThisType();

		static foreach(i; 0..dim) {
			mixin("ret[i][] = data[i][]" ~ op ~ "rhs[i][];");
		}

		return ret;
	}

	@trusted @nogc auto opBinary(string op)(ThisType rhs) immutable pure nothrow {
		auto ret = ThisType();

		static foreach(i; 0..dim) {
			mixin("ret[i][] = data[i][]" ~ op ~ "rhs[i][];");
		}
		
		return ret;
	}

	@trusted @nogc immutable(ThisType) opBinary(string op)(immutable ThisType rhs) immutable pure nothrow {
		auto ret = ThisType();

		static foreach(i; 0..dim) {
			mixin("ret[i][] = data[i][]" ~ op ~ "rhs[i][];");
		}
		
		return ret;
	}

	@trusted @nogc immutable(ThisType) opBinary(string op)(immutable T[size] rhs) immutable pure nothrow {
		auto ret = ThisType();

		static foreach(i; 0..dim) {
			mixin("ret[i][] = data[i][]" ~ op ~ "rhs[];");
		}
		
		return ret;
	}

	@trusted @nogc immutable(ThisType) opBinary(string op)(immutable ForeachType!T rhs) immutable pure nothrow {
		auto ret = ThisType();

		static foreach(i; 0..dim) {
			mixin("ret[i][] = data[i][]" ~ op ~ "rhs;");
		}
		
		return ret;
	}

	@trusted @nogc auto opBinary(string op)(T[size] rhs) pure nothrow {
		auto ret = ThisType();

		static foreach(i; 0..dim) {
			mixin("ret[i][] = data[i][]" ~ op ~ "rhs[];");
		}

		return ret;
	}

	@trusted @nogc auto opBinary(string op)(ref Vec rhs) immutable pure nothrow {
		auto ret = ThisType();

		static foreach(i; 0..dim) {
			mixin("ret[i][] = data[i][]" ~ op ~ "rhs[i];");
		}
		
		return ret;
	}

	@trusted @nogc auto opBinary(string op)(Vec rhs) pure nothrow {
		auto ret = ThisType();

		static foreach(i; 0..dim) {
			mixin("ret[i][] = data[i][]" ~ op ~ "rhs[i];");
		}

		return ret;
	}

	@trusted @nogc auto opBinary(string op)(immutable Vec rhs) immutable pure nothrow {
		auto ret = ThisType();

		static foreach(i; 0..dim) {
			mixin("ret[i][] = data[i][]" ~ op ~ "rhs[i];");
		}

		return ret;
	}

	@trusted @nogc auto opBinaryRight(string op)(ref Vec lhs) immutable pure nothrow {
		auto ret = ThisType();

		static foreach(i; 0..dim) {
			mixin("ret[i][] = lhs[i]" ~ op ~ "data[i][];");
		}

		return ret;
	}

	@trusted @nogc auto opBinaryRight(string op)(Vec lhs) pure nothrow {
		auto ret = ThisType();

		static foreach(i; 0..dim) {
			mixin("ret[i][] = lhs[i]" ~ op ~ "data[i][];");
		}

		return ret;
	}

	/+@trusted @nogc immutable(ThisType) opBinaryRight(string op, L)(L lhs) immutable pure nothrow if((op == "*" || op == "/") && !isArray!L) {
		auto ret = ThisType();

		static foreach(i; 0..dim) {
			mixin("ret[i][] = lhs" ~ op ~ "data[i][];");
		}

		return ret;
	}

	@trusted @nogc ref auto opBinaryRight(string op, T)(T lhs) pure nothrow if((op == "*" || op == "/") && !isArray!T) {
		auto ret = ThisType();

		static foreach(i; 0..dim) {
			mixin("ret[i][] = lhs" ~ op ~ "data[i][];");
		}
		
		return ret;
	}+/

	@trusted @nogc immutable(ThisType) opBinaryRight(string op, size_t size)(immutable T lhs) immutable pure nothrow {
		auto ret = ThisType();

		static foreach(i; 0..dim) {
			mixin("ret[i][] = lhs[]" ~ op ~ "data[i][];");
		}

		return ret;
	}

	@trusted @nogc ref auto opBinaryRight(string op, size_t size)(auto ref T lhs) pure nothrow {
		auto ret = ThisType();

		static foreach(i; 0..dim) {
			mixin("ret[i][] = lhs[]" ~ op ~ "data[i][];");
		}

		return ret;
	}

	@trusted @nogc auto opBinary(string op)(ref Vector!(dim, sizediff_t) rhs) pure nothrow {
		auto ret = ThisType();

		static foreach(i; 0..dim) {
			mixin("ret[i][] = data[i][]" ~ op ~ "rhs[i];");
		}
		
		return ret;
	}

	@trusted @nogc auto opBinary(string op)(Vector!(dim, sizediff_t) rhs) pure nothrow {
		auto ret = ThisType();

		static foreach(i; 0..dim) {
			mixin("ret[i][] = data[i][]" ~ op ~ "rhs[i];");
		}

		return ret;
	}

	@trusted @nogc auto opAssign(VecArrayT!(dim, T) rhs)
	{
		static foreach(i; 0..dim) {
			data[i][] = rhs[i][];
		}

		return this;
	}

	@trusted @nogc auto magnitude_squared() pure nothrow {
		T mag_2 = data[0][]*data[0][] + data[1][]*data[1][] + data[2][]*data[2][];
		return mag_2;
	}

	@trusted @nogc auto magnitude() immutable pure nothrow {
		//version(LDC) pragma(inline, true);
		T mag = data[0][]*data[0][] + data[1][]*data[1][] + data[2][]*data[2][];

		mag = sqrt(mag);
		return mag;
	}

	@trusted @nogc auto normalize() immutable pure nothrow {
		//version(LDC) pragma(inline, true);
		T mag = data[0][]*data[0][] + data[1][]*data[1][] + data[2][]*data[2][];
		mag = 1./sqrt(mag)[];

		auto ret = ThisType();
		ret[0][] = data[0][]*mag[];
		ret[1][] = data[1][]*mag[];
		ret[2][] = data[2][]*mag[];

		return ret;
	}

	@trusted @nogc auto dot(immutable Vec rhs) immutable pure nothrow {
		T ret;
		ret[] = data[0][]*rhs[0] + data[1][]*rhs[1] + data[2][]*rhs[2];
		return ret;
	}

	@trusted @nogc auto dot(ref Vec rhs) immutable pure nothrow {
		T ret;
		ret[] = data[0][]*rhs[0] + data[1][]*rhs[1] + data[2][]*rhs[2];
		return ret;
	}

	@trusted @nogc auto dot(ref Vec rhs) pure nothrow {
		T ret;
		ret[] = data[0][]*rhs[0] + data[1][]*rhs[1] + data[2][]*rhs[2];
		return ret;
	}

	@trusted @nogc auto dot(immutable Vec rhs) pure nothrow {
		T ret;
		ret[] = data[0][]*rhs[0] + data[1][]*rhs[1] + data[2][]*rhs[2];
		return ret;
	}

	/+
	@trusted @nogc auto dot(VecIn)(immutable VecIn rhs) immutable pure nothrow {
		T ret;
		ret[] = data[0][]*rhs[0] + data[1][]*rhs[1] + data[2][]*rhs[2];
		return ret;
	}
	+/
	/+
	/++
	+/
	@trusted @nogc T[size] dot(alias vector)() immutable pure nothrow {
		//version(LDC) pragma(inline, true);
		T[size] val = data[0][] * vector[0] + data[1][] * vector[1] + data[2][] * vector[2];
		return val;
	}

	/++
	+/
	@trusted @nogc T[size] dot(alias vector)() pure nothrow {
		//version(LDC) pragma(inline, true);
		T[size] val = data[0][] * vector[0] + data[1][] * vector[1] + data[2][] * vector[2];
		return val;
	}
	+/
}

/+
auto sqrt(T)(auto ref T data) if(isStaticArray!T) {
	version(LDC) pragma(inline, true);
	static import std.math;
	import std.array : staticArray;

	T ret_data =
		data[]
			.map!(a => std.math.sqrt(a))
			.staticArray!(T.length);

	return ret_data;
}

auto abs(T)(auto ref T data) if(isStaticArray!T) {
	version(LDC) pragma(inline, true);
	static import std.math;
	import std.array : staticArray;

	T ret_data =
		data[]
			.map!(a => std.math.abs(a))
			.staticArray!(T.length);

	return ret_data;
}


auto acos(T)(auto ref T val) {
	//version(LDC) pragma(inline, true);
	static import std.math;
	static if(isStaticArray!T) {
		T res;
		
		foreach(i; 0..T.length) {
			res[i] = std.math.acos(val[i]);
		}
		return res;
	} else {
		return std.math.acos(val);
	}
}
+/
