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
	import std.conv : to;
	import std.stdio : writeln;
	writeln("hello int");

	import opencopter.airfoilmodels;
	
    class DummyAF: AirfoilModel {
        string name;
        this(string _name) {
            name = _name;
        }

        override Chunk get_Cl(Chunk alpha_query, Chunk mach_query) {
            Chunk c;
            return c;
        }

        override Chunk get_Cd(Chunk alpha_query, Chunk mach_query) {
            Chunk c;
            return c;
        }

        override double get_Cl(double alpha_query, double mach_query) {
            return 0.0;
        }

        override double get_Cd(double alpha_query, double mach_query) {
            return 0.0;
        }
    }

    AirfoilModel af1 = new DummyAF("af1");
    AirfoilModel af2 = new DummyAF("af2");
    AirfoilModel af3 = new DummyAF("af3");
    AirfoilModel af4 = new DummyAF("af4");

    auto af_models = [af1, af2, af3, af4];

    size_t[2][] extents = [
        [0, 6],
        [7, 22],
        [23, 26],
        [27,31]
    ];

    auto blade_af = new BladeAirfoil(af_models, extents);

    writeln("blade_af.airfoil_models.length: ", blade_af.airfoil_models.length);
    foreach(blade_af_models; blade_af.airfoil_models) {
        writeln("blade_af_models.length: ", blade_af_models.length);
        foreach(af_model; blade_af_models) {
            writeln("af_model.name: ", af_model.to!DummyAF.name);
        }
        writeln;
    }

	Chunk generate_radius_points(size_t n_sections) {
		import std.conv : to;
		import numd.utility : linspace;
    	auto i = linspace(0., (n_sections-1).to!double, n_sections);
    	Chunk r_array = 1./n_sections.to!double*(i[] + 0.5)[];
    	return r_array;
	}

	double AR = 15;
	auto blade = BladeGeometry(8, 0, 1/AR, blade_af);
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

@nogc Chunk steady_lift_model(immutable Chunk u_p, immutable Chunk u_t, immutable Chunk C_l, immutable Chunk chord) {

	version(LDC) pragma(inline, true);

	Chunk lift_coefficient;

	import std.math : PI;

	immutable Chunk sigma_hat = 2.0*chord[]/(PI*PI);

	immutable Chunk u_squared = u_p[]*u_p[] + u_t[]*u_t[];

	lift_coefficient[] = 0.5*sigma_hat[]*u_squared[]*C_l[];

	return lift_coefficient;
}
