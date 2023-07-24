module opencopter.math.complex;

import opencopter.aircraft;
import opencopter.config;
import opencopter.math;
import opencopter.memory;

import std.algorithm;
import std.array;
import std.complex;
import std.math;

immutable I = complex(0, 1);
immutable Iv = ichunk(0, 1);
immutable Z = zero;
immutable one = ichunk(1, 0);

@nogc Chunk re(ref IChunk x) {
	return x[].map!(a => a.re).staticArray!Chunk;
}

@nogc Chunk re(IChunk x) {
	return x[].map!(a => a.re).staticArray!Chunk;
}

@nogc Chunk im(ref IChunk x) {
	return x[].map!(a => a.im).staticArray!Chunk;
}

@nogc Chunk im(IChunk x) {
	return x[].map!(a => a.im).staticArray!Chunk;
}

@nogc Chunk zero() {
	Chunk r;
	r[] = 0;
	return r;
}

@nogc IChunk ichunk(Chunk r, double i) {
	IChunk res;
	foreach(idx; 0..Chunk.length) {
		res[idx].re = r[idx];
		res[idx].im = i;
	}
	return res;
}

@nogc IChunk ichunk(double r, Chunk i) {
	IChunk res;
	foreach(idx; 0..Chunk.length) {
		res[idx].re = r;
		res[idx].im = i[idx];
	}
	return res;
}

@nogc IChunk ichunk(double r, double i) {
	IChunk res;
	foreach(idx; 0..Chunk.length) {
		res[idx].re = r;
		res[idx].im = i;
	}
	return res;
}

@nogc IChunk ichunk(Chunk r, Chunk i) {
	IChunk res;
	foreach(idx; 0..Chunk.length) {
		res[idx].re = r[idx];
		res[idx].im = i[idx];
	}
	return res;
}

@nogc IChunk log(IChunk x) {
	static import std.complex;
	IChunk res;
	foreach(i; 0..IChunk.length) {
		res[i] = std.complex.log(x[i]);
	}
	return res;
}

@nogc IChunk log10(IChunk x) {
	static import std.complex;
	IChunk res;
	foreach(i; 0..IChunk.length) {
		res[i] = std.complex.log10(x[i]);
	}
	return res;
}

@nogc IChunk pow(IChunk x, double n) {
	static import std.complex;
	IChunk res;
	foreach(i; 0..IChunk.length) {
		res[i] = std.complex.pow(x[i], n);
	}
	return res;
}

@nogc IChunk sqrt(IChunk x) {
	static import std.complex;
	IChunk res;
	foreach(i; 0..IChunk.length) {
		res[i] = std.complex.sqrt(x[i]);
	}
	return res;
}

/+@nogc IChunk atan(IChunk x) {
	static import std.complex;
	IChunk res;
	foreach(i; 0..IChunk.length) {
		res[i] = std.complex.atan(x[i]);
	}
	return res;
}+/
