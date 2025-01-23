module opencopter.airfoilmodels;

public import opencopter.airfoilmodels.aerodas;
public import opencopter.airfoilmodels.c81;
public import opencopter.airfoilmodels.thinaf;

import opencopter.config;
import opencopter.memory;

import std.exception;

class AirfoilModel {
    Chunk get_Cl(Chunk alpha_query, Chunk mach_query) { assert(0); }
    Chunk get_Cd(Chunk alpha_query, Chunk mach_query) { assert(0); }
    double get_Cl(double alpha_query, double mach_query) { assert(0); }
    double get_Cd(double alpha_query, double mach_query) { assert(0); }
    double lift_curve_slope() { assert(0); }
    double zero_lift_aoa() { assert(0); }
    //Chunk get_Cl(Chunk alpha_query, Chunk mach_query);
    //Chunk get_Cd(Chunk alpha_query, Chunk mach_query);
    //double get_Cl(double alpha_query, double mach_query);
    //double get_Cd(double alpha_query, double mach_query);
}

struct AirfoilState {
    Chunk C_l;
    Chunk C_d;
}

unittest {
    import std.conv : to;
    import std.stdio : writeln;

    writeln("hello");
    class DummyAF: AirfoilModel {
        string name;
        this(string _name) {
            name = _name;
        }

        override double lift_curve_slope() {
            return 2.0*PI;
        }

        override double zero_lift_aoa() {
            return 0;
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
}

class BladeAirfoil {
    alias AirfoilModels = AirfoilModel[];
    AirfoilModels[] airfoil_models;

    this(AirfoilModel[] af_models, size_t[2][] extents) {
        import std.stdio : writeln;
        enforce(af_models.length == extents.length, "The supplied number of airfoil models should be the same number of supplied extents.");

        immutable size_t num_chunks = (extents[$-1][1] - extents[0][0] + 1)/chunk_size;

        airfoil_models = new AirfoilModels[num_chunks];

        debug writeln("num_chunks: ", num_chunks);

        foreach(e_idx, extent; extents) {
            debug writeln("extent: ", extent);
            foreach(c_idx; (e_idx != 0 ? extent[0] - 1 : extent[0])..(e_idx != extents.length - 1 ? extent[1] : extent[1] + 1)) {
                debug writeln("c_idx: ", c_idx);
                immutable chunk_idx = c_idx/chunk_size;

                debug writeln("chunk_idx: ", chunk_idx);

                airfoil_models[chunk_idx] ~= af_models[e_idx];

            }
        }
    }

    AirfoilState compute_coeffiecients(size_t chunk_idx, Chunk aoa, Chunk mach) {
        AirfoilState state;
        auto chunk_models = airfoil_models[chunk_idx];
        if(chunk_models.length == 1) {
            state.C_l[] = chunk_models[0].get_Cl(aoa, mach)[];
            state.C_d[] = chunk_models[0].get_Cd(aoa, mach)[];
        } else {
            foreach(c_idx; 0..chunk_size) {
                state.C_l[c_idx] = chunk_models[c_idx].get_Cl(aoa[c_idx], mach[c_idx]);
                state.C_d[c_idx] = chunk_models[c_idx].get_Cd(aoa[c_idx], mach[c_idx]);
            }
        }
        return state;
    }

    Chunk lift_curve_slope(size_t chunk_idx) {
        Chunk Cl_alpha;
        auto chunk_models = airfoil_models[chunk_idx];
        if(chunk_models.length == 1) {
            Cl_alpha[] = chunk_models[0].lift_curve_slope();
        } else {
            foreach(c_idx; 0..chunk_size) {
                Cl_alpha[c_idx] = chunk_models[c_idx].lift_curve_slope();
            }
        }

        return Cl_alpha;
    }

    Chunk zero_lift_aoa(size_t chunk_idx) {
        Chunk Cl_alpha;
        auto chunk_models = airfoil_models[chunk_idx];
        if(chunk_models.length == 1) {
            Cl_alpha[] = chunk_models[0].zero_lift_aoa();
        } else {
            foreach(c_idx; 0..chunk_size) {
                Cl_alpha[c_idx] = chunk_models[c_idx].zero_lift_aoa();
            }
        }

        return Cl_alpha;
    }
}
