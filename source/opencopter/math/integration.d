module opencopter.math.integration;

import opencopter.aircraft;
import opencopter.config;
import opencopter.math;
import opencopter.memory;

import std.math;

double integrate_trapaziodal(string value, BS, BG)(auto ref BS blade, auto ref BG blade_geom) {
	double val = 0;
	import std.stdio : writeln;

	foreach(size_t idx, ref BladeStateChunk chunk; blade.chunks) {
		Chunk dr;
		Chunk tmp_sum;
		mixin("tmp_sum[0..$-1] = 0.5*(chunk."~value~"[1..$] + chunk."~value~"[0..$-1]);");

		dr[0..$-1] = blade_geom.chunks[idx].r[1..$] - blade_geom.chunks[idx].r[0..$ - 1];
		if(idx != blade.chunks.length - 1) {
			dr[$-1] = blade_geom.chunks[idx + 1].r[0] - blade_geom.chunks[idx].r[$-1];
			mixin("tmp_sum[$-1] = 0.5*(blade.chunks[idx + 1]."~value~"[0] + chunk."~value~"[$-1]);");
		} else {
			tmp_sum[$-1] = 0;
			dr[$-1] = 0;
		}

		tmp_sum[] *= abs(dr)[];

		import std.algorithm : sum;
		val += tmp_sum[].sum;
	}

	return val;
}

double integrate_trapaziodal(BG)(ref Chunk[] f, auto ref BG blade_geom) {
	double val = 0;
	import std.stdio : writeln;

	foreach(size_t idx, ref Chunk chunk; f) {
		Chunk dr;
		Chunk tmp_sum;
		tmp_sum[0..$-1] = 0.5*(chunk[1..$] + chunk[0..$-1]);

		dr[0..$-1] = blade_geom.chunks[idx].r[1..$] - blade_geom.chunks[idx].r[0..$ - 1];
		if(idx != f.length - 1) {
			dr[$-1] = blade_geom.chunks[idx + 1].r[0] - blade_geom.chunks[idx].r[$-1];
			tmp_sum[$-1] = 0.5*(f[idx + 1][0] + chunk[$-1]);
		} else {
			tmp_sum[$-1] = 0;
			dr[$-1] = 0;
		}

		tmp_sum[] *= abs(dr)[];

		import std.algorithm : sum;
		val += tmp_sum[].sum;
	}

	return val;
}

double integrate_trapaziodal(ref Chunk[] f, ref Chunk[] r) {
	double val = 0;
	import std.stdio : writeln;

	foreach(size_t idx, ref Chunk chunk; f) {
		Chunk dr;
		Chunk tmp_sum;
		tmp_sum[0..$-1] = 0.5*(chunk[1..$] + chunk[0..$-1]);

		dr[0..$-1] = r[idx][1..$] - r[idx][0..$ - 1];
		if(idx != f.length - 1) {
			dr[$-1] = r[idx + 1][0] - r[idx][$-1];
			tmp_sum[$-1] = 0.5*(f[idx + 1][0] + chunk[$-1]);
		} else {
			tmp_sum[$-1] = 0;
			dr[$-1] = 0;
		}

		tmp_sum[] *= abs(dr)[];

		import std.algorithm : sum;
		val += tmp_sum[].sum;
	}

	return val;
}

double integrate_trapaziodal_reverse(ref Chunk[] f, ref Chunk[] r) {
	double val = 0;
	import std.stdio : writeln;

	foreach(size_t idx, ref Chunk chunk; f) {
		Chunk dr;
		Chunk tmp_sum;
		tmp_sum[0..$-1] = 0.5*(chunk[1..$] + chunk[0..$-1]);

		dr[0..$-1] = r[idx][0..$ - 1] - r[idx][1..$];
		if(idx != f.length - 1) {
			dr[$-1] = r[idx][$-1] - r[idx + 1][0];
			tmp_sum[$-1] = 0.5*(f[idx + 1][0] + chunk[$-1]);
		} else {
			tmp_sum[$-1] = 0;
			dr[$-1] = 0;
		}

		tmp_sum[] *= abs(dr)[];

		import std.algorithm : sum;
		val += tmp_sum[].sum;
	}

	return val;
}

double integrate_trapaziodal(double[] f, double[] r) {
	double sum = 0;
	import std.math : abs;
	//import std.range : iota;
	//import std.algorithm;
	foreach(size_t idx; 1..f.length) {
		sum += (f[idx - 1] + f[idx])*abs(r[idx] - r[idx - 1]);
	}
	//return 0.5*iota(1, f.length).map!(idx => (f[idx - 1] + f[idx])*(r[idx] - r[idx - 1])).sum;

	return 0.5*sum;
}

unittest {
	import numd.utility;

	import std.algorithm;
	import std.array;
	import std.conv;
	import std.math;
	import std.stdio;
	import std.range;
	
	writeln("hello int");
	void set_array(ref Chunk[] chunk_data, double[] data) {

		foreach(c_idx, ref chunk; chunk_data) {
			immutable out_start_idx = c_idx*chunk_size;

			immutable remaining = data.length - out_start_idx;
			
			immutable out_end_idx = remaining > chunk_size ? (c_idx + 1)*chunk_size : out_start_idx + remaining;
			immutable in_end_idx = remaining > chunk_size ? chunk_size : remaining;

			chunk[0..in_end_idx] = data[out_start_idx..out_end_idx];
		}
	}

	double[] generate_radius_points(double end, size_t n_sections) {
	    return iota(1.0, n_sections + 1.0).map!((n) {
	    	immutable psi = n*PI/(n_sections.to!double + 1.0);
	    	auto r = end*0.5*(cos(psi) + 1.0).to!double;
	    	return r;
	    }).retro.array;
	}

	

	size_t len = 512;
	auto _r = generate_radius_points(PI, len);
	writeln(_r);
	double[] _f = new double[len];

	foreach(i; 0..len) {
		_f[i] = sin(_r[i]);
	}

	Chunk[] r = new Chunk[len/chunk_size];
	Chunk[] f = new Chunk[len/chunk_size];
	set_array(r, _r);
	set_array(f, _f);

	auto val = integrate_trapaziodal(f, r);
	writeln("integrated val: ", val);

	val = integrate_trapaziodal(_f, _r);
	writeln("integrated val: ", val);
}
