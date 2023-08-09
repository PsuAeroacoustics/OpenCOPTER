module opencopter.airfoilmodels.thinaf;

import opencopter.airfoilmodels;
import opencopter.memory;


import std.math;

class ThinAirfoil: AirfoilModel {

	double C_l_alpha_0;

	this(double _C_l_alpha_0) {
		C_l_alpha_0 = _C_l_alpha_0;
	}

    override Chunk get_Cl(Chunk alpha_query, Chunk mach_query) {
		Chunk C_l = 2.0*PI*alpha_query[] + C_l_alpha_0;
		return C_l;
	}

    override Chunk get_Cd(Chunk alpha_query, Chunk mach_query) {
		Chunk C_d = 0.0;
		return C_d;
	}

    override double get_Cl(double alpha_query, double mach_query) {
		return 2.0*PI*alpha_query + C_l_alpha_0;
	}

    override double get_Cd(double alpha_query, double mach_query) {
		return 0.0;
	}
}
