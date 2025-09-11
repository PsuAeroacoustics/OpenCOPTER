module opencopter.memory;

import opencopter.config;

import std.complex : Complex;

alias Chunk = double[chunk_size];
alias IChunk = Complex!(double)[chunk_size];

double[][] allocate_dense(size_t r, size_t c) {
    double[] _f = new double[r*c];
    double[][] f;
    f.length = r;
    foreach(idx; 0..r) {
        f[idx] = _f[idx*c..idx*c+c];
    }
	return f;
}

Chunk[][] allocate_dense_chunk(size_t r, size_t c) {
    Chunk[] _f = new Chunk[r*c];
    Chunk[][] f;
    f.length = r;
    foreach(idx; 0..r) {
        f[idx] = _f[idx*c..idx*c+c];
    }
	return f;
}

Chunk[][] allocate_dense_chunk_aliased(size_t r, size_t c) {
	immutable col_chunks = c/chunk_size;
    Chunk[] _f = new Chunk[r*col_chunks];
    Chunk[][] f;
    f.length = r;
    foreach(idx; 0..r) {
        f[idx] = _f[idx*col_chunks..idx*col_chunks+col_chunks];
    }
	return f;
}

/++
 +	This exists purely to dirty a function signature
 +	so it can be used with the python wrapper. If a function
 +	calls this it will remove pure and @system that is auto
 +	applied to functions.
 +/
@system extern (C) void nop() {

}


/++
 +	Array Container that exists purely to interop with python.
 +	Interfacing arrays with python don't behave the same way as
 +	D. We end up needing to return pointers to the array elements
 +	to avoid the constant construction and destruction of stored
 +	structs. 
 +/
struct Array(__T) {
	alias T = __T;
	T[] data;

	size_t curr_idx;

	this(size_t size) {
		data = new T[size];
		curr_idx = 0;
	}

	this(T[] slice) {
		data = slice;
		curr_idx = 0;
	}

	version(Have_pyd) {
		import pyd.pyd : d_to_python;

		typeof(this)* __iter__() {
			nop;
			curr_idx = 0;
			return &this;
		}

		import deimos.python.Python : PyObject;
		PyObject* next() {
			nop;
			if(curr_idx < data.length) {
				curr_idx++;
				return d_to_python(&(data[curr_idx - 1]));
			} else {
				return null;
			}
		}

		size_t __len__() {
			nop;
			return data.length;
		}

		size_t len() {
			nop;
			return data.length;
		}
	}

	string toString() {
		import std.format : format;
		return format!"%s"(data);
	}

	bool empty() {
		return curr_idx >= data.length;
	}

	void popFront() {
		curr_idx++;
	}

	T* front() {
		return &data[curr_idx];
	}

	auto opOpAssign(string op, T)(T value)
	{
		static if(op == "~") {
			data ~= value;
		}
		return this;
	}

	auto opOpAssign(string op, T)(T* value)
	{
		static if(op == "~") {
			data ~= *value;
		}
		return this;
	}

	auto opAssign(T[] slice) {
		data = slice;
		return this;
	}

	auto opAssign(ref Array!T arr) {
		data = arr.data;
		return this;
	}

	auto opAssign(Array!T arr) {
		data = arr.data;
		return this;
	}

	import std.traits : isBasicType;
	static if(isBasicType!T) {
		pragma(msg, "is basic: ", T.stringof);
		T opIndex(size_t i) {
			nop;
			assert(i < data.length);
			return data[i];
		}
	} else {
		pragma(msg, "is not basic: ", T.stringof);
		T* opIndex(size_t i) {
			nop;
			assert(i < data.length);
			return &(data[i]);
		}
	}

 
	void opIndexAssign(T* t, size_t i) {
		assert(i < data.length);
		data[i] = *t;
	}

	void opIndexAssign(T t, size_t i) {
		assert(i < data.length);
		data[i] = t;
	}

	size_t length() {
		nop;
		return data.length;
	}

	int opApply(int delegate(ref T) dg) {
		foreach(i; 0..length) {
			dg(data[i]);
		}
		return 1;
	}

	int opApply(int delegate(ref size_t, ref T) dg) {
		foreach(i; 0..length) {
			dg(i, data[i]);
		}
		return 1;
	}

	auto ref opSlice(size_t i, size_t j) {
		return data[i..j];
	}

	auto ref opSlice() {
		return data;
	}

	size_t opDollar() {
		return data.length;
	}
}

auto ptr_at_idx(A)(auto ref A array, size_t i) {
	import std.traits : isInstanceOf;
	pragma(msg, A);
	static if(isInstanceOf!(Array, A)) {
		pragma(msg, "is Array");
		return array[i];
	} else {
		pragma(msg, "is D array");
		return &array[i];
	}
}

auto at_idx(A)(auto ref A array, size_t i) {
	import std.traits : isInstanceOf;
	pragma(msg, A);
	static if(isInstanceOf!(Array, A)) {
		pragma(msg, "is Array");
		return *array[i];
	} else {
		pragma(msg, "is D array");
		return array[i];
	}
}

enum ArrayContainer : int {
	none,
	array
}

template ArrayDeclMixin(ArrayContainer AC, Type, string name) {
	static if(AC == ArrayContainer.none) {
		// Nitya: String mixin enable string constants to be compiled
		//        as regular D code and inserted into the program
		mixin("Type[] "~name~";");
	} else {
		//import opencopter.python : Array;
		mixin("Array!(Type) "~name~";");
	}
}

string array_ctor_mixin(ArrayContainer AC, string Type, string name, string size)() {
	static if(AC == ArrayContainer.none) {
		return name~" = new "~Type~"["~size~"];";
	} else {
		return "/+import opencopter.memory : Array;+/\n"~name~" = Array!("~Type~")("~size~");";
	}
}

string array_ctor_mixin_slice(ArrayContainer AC, string Type, string name, string slice)() {
	static if(AC == ArrayContainer.none) {
		return name~" = "~slice~";";
	} else {
		return "/+import opencopter.memory : Array;+/\n"~name~" = Array!("~Type~")("~slice~");";
	}
}

void set_array(ref Chunk[] chunks, double[] data) {
	foreach(c_idx, ref chunk; chunks) {
		immutable out_start_idx = c_idx*chunk_size;

		immutable remaining = data.length - out_start_idx;
		// Nitya: RESULT = (COND) ? ( STATEMENT IF TRUE) : (STATEMENT IF FALSE)
		immutable out_end_idx = remaining > chunk_size ? (c_idx + 1)*chunk_size : out_start_idx + remaining;
		immutable in_end_idx = remaining > chunk_size ? chunk_size : remaining;

		chunk[0..in_end_idx] = data[out_start_idx..out_end_idx];
	}
}

void set_array(ref Array!Chunk chunks, double[] data) {
	foreach(c_idx, ref chunk; chunks) {
		immutable out_start_idx = c_idx*chunk_size;

		immutable remaining = data.length - out_start_idx;
		
		immutable out_end_idx = remaining > chunk_size ? (c_idx + 1)*chunk_size : out_start_idx + remaining;
		immutable in_end_idx = remaining > chunk_size ? chunk_size : remaining;

		chunk[0..in_end_idx] = data[out_start_idx..out_end_idx];
	}
}
