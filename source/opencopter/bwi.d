module opencopter.bwi;

import opencopter.aircraft;
import opencopter.config;
import opencopter.memory;
import opencopter.math;
import opencopter.wake;
import opencopter.bladeelement;
import opencopter.aircraft;

import std.math;

immutable size_t nPoints = 4;

extern (C++) struct BWIinputsChunk {
    Chunk miss_dist;
    Chunk r_c_ave;
    Chunk gamma_w;
	Chunk gamma_sec;
    Chunk theta_v;
    Chunk phi_v;
    // Angle of interaction!! 
    //double wake_age;
    //double azimuth;
    //Vec3 bwi_points;
    // Nitya: Do I define a constructor?? - We don't!
    // InteractionPoints interaction_pts;
}

extern (C++) struct InteractionPoints {

    double[nPoints] x_pos;
    double[nPoints] y_pos;
    double[nPoints] z_pos;
    double[nPoints] miss_dist;
    double[nPoints] r_c;
    double[nPoints] gamma_w;
    double[nPoints] gamma_sec;

}

extern (C++) struct TipVortexInteractionT(ArrayContainer AC){
    import std.stdio : writeln;
    mixin ArrayDeclMixin!(AC, BWIinputsChunk, "BWI_inputs");

    size_t length;

	this(size_t wake_history) {

		assert(wake_history % chunk_size == 0);
		immutable num_chunks = wake_history/chunk_size;
        
		mixin(array_ctor_mixin!(AC, "BWIinputsChunk", "BWI_inputs", "num_chunks"));
        debug writeln("BWI_inputs:", BWI_inputs.length);

        foreach(ref chunk; BWI_inputs) {
		 	chunk.gamma_w[] = 0;
		    chunk.r_c_ave = 0.00001;
		}
        length = 0;
	}
    
}

alias TipVortexInteraction = TipVortexInteractionT!(ArrayContainer.none);

double[] get_BWIinputs(string value, BWI)(auto ref BWI VortexInteraction) {
	
	immutable elements = VortexInteraction.BWI_inputs.length * chunk_size;

	double[] d = new double[elements];

	foreach(c_idx, ref chunk; VortexInteraction.BWI_inputs) {
		immutable out_start_idx = c_idx*chunk_size;

		immutable remaining = elements - out_start_idx;
		
		immutable out_end_idx = remaining > chunk_size ? (c_idx + 1)*chunk_size : out_start_idx + remaining;
		immutable in_end_idx = remaining > chunk_size ? chunk_size : remaining;

		mixin("d[out_start_idx..out_end_idx] = chunk."~value~"[0..in_end_idx];");
	}
	return d;
}

void calculate_BWI_points (BWI_inputs, W, BS)(auto ref BWI_inputs bwi_inputs, auto ref W wake, auto ref BS blade_state, size_t rotor_idx, size_t blade_idx){
    // What we need
    // 1. Wake: x, y, z 
    // 2. r_c and miss distance
    // 3. local airfoil points: vehicle.aircraft.rotors.blades.chuncks (vtk.d ln340), it's calculated in bladeelement.d too!
    // Wake history is stored every 1 deg., therefore 1st interaction should be around 90 deg. 

    double[] x_vortex;
    double[] y_vortex;
    //double[] z_vortex;
    double[] x_blade;
    double[] y_blade;
    //double[] z_blade;

    immutable double epsilon = 0.0001;
    auto current_tip_filament = wake.rotor_wakes[rotor_idx].tip_vortices[blade_idx];

    x_blade = get_wake_component!"x"(blade_state);
    y_blade = get_wake_component!"y"(blade_state);
    //z_blade = get_state_array!z(blade_state);

    x_vortex = get_wake_component!"x"(current_tip_filament); 
    y_vortex = get_wake_component!"y"(current_tip_filament);

    size_t idx = 0;
    size_t len;
    len = x_blade.length;
    /* foreach (i_point; 0..x_vortex.length){
        //while (idx < nPoints) {
        if(x_vortex[i_point].isNaN) {
            break;
        }
            if((abs(x_blade[$]-x_vortex[i_point])<epsilon) && (abs(y_blade[$]-y_vortex[i_point])<epsilon)) {
                idx++;
            }
        
        //}
    }  */
    

    foreach (wake_idx; 0..x_vortex.length){
        if(x_vortex[wake_idx].isNaN) {
            break;
        }
        foreach(span_idx; 0..x_blade.length){
            if((abs(x_blade[span_idx]-x_vortex[wake_idx])<epsilon) && (abs(y_blade[span_idx]-y_vortex[wake_idx])<epsilon)) {
                idx++;
                break;
            }
        }
        if(idx >= nPoints){
            break;
        }
    }

} 