module opencopter.liftmodels.steadylift;

import opencopter.aircraft;
import opencopter.atmosphere;
import opencopter.inflow;
import opencopter.math;
import opencopter.memory;
import opencopter.trim;

import std.algorithm;
import std.array;
import std.math;

unittest {
	Chunk generate_radius_points(size_t n_sections) {
		import std.conv : to;
		import numd.utility : linspace;
    	auto i = linspace(0., (n_sections-1).to!double, n_sections);
    	Chunk r_array = 1./n_sections.to!double*(i[] + 0.5)[];
    	return r_array;
	}

	double AR = 15;
	auto blade = BladeGeometry(8, 0, 1/AR);
	//BladeGeometryChunk chunk;
	blade.chunks[0].r[] = generate_radius_points(8);// linspace(0., 1., chunk_size);
	blade.chunks[0].twist[] = (12*(PI/180.0))/blade.chunks[0].r[];
	blade.chunks[0].chord[] = 1/AR;
	blade.chunks[0].sweep[] = 0;
	blade.chunks[0].C_l_alpha[] = 2.0*PI;
	blade.chunks[0].alpha_0[] = 0;
	double rotor_radius = 4;


	double C_T = 0.008;
	double lambda = sqrt(0.5*C_T);
	Chunk inflow = lambda;
	Chunk u_t = blade.chunks[0].r[];

	Chunk inflow_angle = atan2(inflow, u_t)[];

	auto bs = BladeState(1, blade);

	import std.stdio : writeln;
	writeln;

	bs.chunks[0].dC_T = 0;//steady_lift_model(inflow, rotor_radius, blade.chunks[0], u_t, inflow_angle, 0);

	immutable blade_C_T = integrate_trapaziodal!"dC_T"(bs, blade);

	writeln("C_T: ", blade_C_T);
	writeln("u_t: ", u_t);
	writeln("inflow_angle: ", inflow_angle);
	writeln("inflow: ", inflow);
	writeln("lift_coefficient: ", bs.chunks[0].dC_T);
	writeln;
}

@nogc Chunk steady_lift_model(BC, BSC)(Chunk inflow, double rotor_radius, auto ref BC blade_chunk, auto ref BSC blade_state_chunk, immutable Chunk u_t, immutable Chunk inflow_angle, double time) if(is_blade_geometry_chunk!BC) {

	version(LDC) pragma(inline, true);

	Chunk lift_coefficient;

	import std.math : PI;
	import std.stdio : writeln;

	//immutable Chunk sigma_hat = 2.0*blade_chunk.chord[]/(PI*PI);
	//immutable Chunk sigma_hat = blade_chunk.chord[]/(PI*rotor_radius);
	immutable Chunk sigma_hat = blade_chunk.chord[]/(2.0*/+PI*+/PI);

	//immutable Chunk corrected_u_t = u_t[].map!(a => a < 0 ? 0 : a).staticArray!Chunk;
	immutable Chunk u_squared = inflow[]*inflow[] + u_t[]*u_t[];
	//immutable Chunk u_squared = (inflow[] + corrected_u_t[])*(inflow[] + corrected_u_t[]);
	//immutable Chunk angle_of_attack = blade_chunk.twist[] - inflow_angle[];
	immutable Chunk C_l = blade_chunk.C_l_alpha[]*blade_state_chunk.aoa[] + blade_chunk.alpha_0[];

	lift_coefficient[] = 0.5*sigma_hat[]*u_squared[]*C_l[];

	return lift_coefficient;
}
