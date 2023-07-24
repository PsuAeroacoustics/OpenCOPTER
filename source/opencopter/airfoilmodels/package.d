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
        [0, 9],
        [10, 22],
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

        size_t last_chunk_idx  = 0;
        size_t last_sub_idx = 0;
        foreach(e_idx, extent; extents) {
            immutable extent_size = extent[1] - extent[0] + 1;

            debug writeln("extent: ", extent);
            debug writeln("extent_size: ", extent_size);
            
            if((extent_size >= chunk_size) && (last_sub_idx == 0)) {
                immutable number_of_chunks = extent_size/chunk_size;
                debug writeln("number_of_chunks: ", number_of_chunks);
                foreach(chunk_idx; 0..number_of_chunks) {
                    airfoil_models ~= [af_models[e_idx]];
                    last_chunk_idx++;
                }

                //debug writeln(airfoil_models);
                immutable remainder = extent_size%chunk_size;

                debug writeln("remainder: ", remainder);
                if(remainder > 0) {
                    airfoil_models.length++;
                    //debug writeln(airfoil_models);
                    foreach(sub_idx; 0..remainder) {
                        airfoil_models[last_chunk_idx] ~= af_models[e_idx];
                        last_sub_idx++;
                    }
                }
            } else if((extent_size >= chunk_size) && (last_sub_idx != 0)) {
                immutable remaining_in_chunk = chunk_size - last_sub_idx;
                debug writeln("remaining in chunk: ", remaining_in_chunk);
                foreach(sub_idx; 0..remaining_in_chunk) {
                    airfoil_models[last_chunk_idx] ~= af_models[e_idx];
                }
                last_sub_idx = 0;
                last_chunk_idx++;

                immutable number_of_chunks = (extent_size - remaining_in_chunk)/chunk_size;

                debug writeln("number_of_chunks: ", number_of_chunks);
                foreach(chunk_idx; 0..number_of_chunks) {
                    airfoil_models ~= [af_models[e_idx]];
                    last_chunk_idx++;
                }

                immutable remainder = (extent_size - remaining_in_chunk)%chunk_size;
                debug writeln("remaining remainder: ", remainder);
                if(remainder > 0) {
                    airfoil_models.length++;
                    foreach(sub_idx; 0..remainder) {
                        airfoil_models[last_chunk_idx] ~= af_models[e_idx];
                        last_sub_idx++;
                    }
                }
            } else {
                foreach(sub_idx; 0..extent_size) {
                    if(last_sub_idx >= chunk_size) {
                        airfoil_models.length++;
                        last_chunk_idx++;
                        last_sub_idx = 0;
                        airfoil_models[last_chunk_idx] ~= af_models[e_idx];
                    } else {
                        airfoil_models[last_chunk_idx] ~= af_models[e_idx];
                    }
                    
                    last_sub_idx++;
                }
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
}
