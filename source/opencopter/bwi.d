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

extern (C++) struct directionVec {
    double x;
    double y;
    double z;
}

extern (C++) struct BWIinputsChunk {
    Chunk miss_dist;
    Chunk r_c_ave;
    Chunk gamma_w;
	Chunk gamma_sec;
    // Chunk theta_v;
    // Chunk phi_v;
    size_t[chunk_size] bladeSec_idx;
    directionVec[chunk_size] r_blade;
    directionVec[chunk_size] r_vortex;

    /* Chunk[] dist;

    this(size_t num_blades){
        dist.length = num_blades;
    } */
    // Angle of interaction!! 
    //double wake_age;
    //double azimuth;
    //Vec3 bwi_points;
    // Nitya: Do I define a constructor?? - We don't!
    // InteractionPoints interaction_pts;
}

extern (C++) struct InteractionPoints {
    size_t[nPoints] wake_idx;
    size_t[nPoints] bladeSec_idx;
    double[nPoints] miss_dist;
    double[nPoints] r_c;
    double[nPoints] gamma_w;
    double[nPoints] gamma_sec;
    directionVec[nPoints] r_blade;
    directionVec[nPoints] r_vortex;
}

extern (C++) struct TipVortexInteractionT(ArrayContainer AC){
    import std.stdio : writeln;
    mixin ArrayDeclMixin!(AC, BWIinputsChunk, "BWI_inputs");

    InteractionPoints interaction_pts;
    size_t length;

	this(size_t wake_history, size_t blade_elements) {

		assert(wake_history % chunk_size == 0);
		immutable num_chunks = wake_history/chunk_size;
        
		mixin(array_ctor_mixin!(AC, "BWIinputsChunk", "BWI_inputs", "num_chunks"));
        debug writeln("BWI_inputs:", BWI_inputs.length);

        length = 0;
	}
    
}

extern (C++) struct VortexInteractionT(ArrayContainer AC){
    mixin ArrayDeclMixin!(AC, TipVortexInteractionT!(AC), "tip_vortex_interaction");

    size_t length;

	this(size_t num_blades, size_t wake_history, size_t blade_elements) {
        
		mixin(array_ctor_mixin!(AC, "TipVortexInteractionT!(AC)", "tip_vortex_interaction", "num_blades"));
        foreach(ref interaction; tip_vortex_interaction) {
			interaction = TipVortexInteractionT!AC(wake_history, blade_elements);
		}
        
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

double[] get_BWIinputs_directionVec(string value, string value2, BWI)(auto ref BWI VortexInteraction){
    immutable elements = VortexInteraction.BWI_inputs.length * chunk_size;

	double[] d = new double[elements];

	foreach(c_idx, ref chunk; VortexInteraction.BWI_inputs) {
		immutable out_start_idx = c_idx*chunk_size;

		immutable remaining = elements - out_start_idx;
		
		immutable out_end_idx = remaining > chunk_size ? (c_idx + 1)*chunk_size : out_start_idx + remaining;
		immutable in_end_idx = remaining > chunk_size ? chunk_size : remaining;

        foreach(idx; 0..in_end_idx){
            mixin("d[out_start_idx + idx] = chunk."~value~"[idx]."~value2~";");
        }
		
	}
	return d;
}

size_t[] get_blade_sec_idx(string value, BWI)(auto ref BWI VortexInteraction) {
	
	immutable elements = VortexInteraction.BWI_inputs.length * chunk_size;

	size_t[] d = new size_t[elements];

	foreach(c_idx, ref chunk; VortexInteraction.BWI_inputs) {
		immutable out_start_idx = c_idx*chunk_size;

		immutable remaining = elements - out_start_idx;
		
		immutable out_end_idx = remaining > chunk_size ? (c_idx + 1)*chunk_size : out_start_idx + remaining;
		immutable in_end_idx = remaining > chunk_size ? chunk_size : remaining;

		mixin("d[out_start_idx..out_end_idx] = chunk."~value~"[0..in_end_idx];");
	}
	return d;
}

size_t[] get_interaction_points(string value, BWI)(auto ref BWI VortexInteraction){
    size_t[] d = new size_t[nPoints];
    mixin("d[0..$] = VortexInteraction.interaction_pts."~value~"[0..$];");
    return d;
}

double[] get_interaction_point_components(string value, BWI)(auto ref BWI VortexInteraction){
    double[] d = new double[nPoints];
    mixin("d[0..$] = VortexInteraction.interaction_pts."~value~"[0..$];");
    return d;
}

double[] get_interaction_point_directionVec(string value, string value2, BWI)(auto ref BWI VortexInteraction){
    double[] d = new double[nPoints];
    foreach(idx;0..nPoints){
        mixin("d[idx] = VortexInteraction.interaction_pts."~value~"[idx]."~value2~";");
    }
    return d;
}

void calculate_BWI_points (W)(auto ref W wake, size_t rotor_idx, size_t blade_idx){
    // What we need
    // 1. Wake: x, y, z 
    // 2. r_c and miss distance
    // 3. local airfoil points: vehicle.aircraft.rotors.blades.chuncks (vtk.d ln340), it's calculated in bladeelement.d too!
    // Wake history is stored every 1 deg., therefore 1st interaction should be around 90 deg. 

    double[] miss_dist;
    double[] r_c;
    double[] gamma_w;
    double[] gamma_sec;
    size_t[] bladeSec_idx;
    size_t i_pt=0;

    foreach (i_blade_idx; 0..wake.rotor_wakes[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction.length){
        miss_dist = get_BWIinputs!"miss_dist"(wake.rotor_wakes[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction[i_blade_idx]);
        r_c = get_BWIinputs!"r_c_ave"(wake.rotor_wakes[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction[i_blade_idx]);
        bladeSec_idx = get_blade_sec_idx!"bladeSec_idx"(wake.rotor_wakes[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction[i_blade_idx]);
        i_pt = 0;
        double local_min = miss_dist[0];
        foreach(idx; 0..miss_dist.length-1){
            if(idx == 0) {
                if(miss_dist[idx]<miss_dist[idx+1]){
                    local_min = miss_dist[idx];
                    if((i_pt < nPoints) && (local_min<4*r_c[idx]*r_c[idx])){
                        wake.rotor_wakes[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction[i_blade_idx].interaction_pts.wake_idx[i_pt] = idx;
                        wake.rotor_wakes[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction[i_blade_idx].interaction_pts.miss_dist[i_pt] = miss_dist[idx];
                        wake.rotor_wakes[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction[i_blade_idx].interaction_pts.r_c[i_pt] = r_c[idx];
                        fill_interactionPts(wake.rotor_wakes[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction[i_blade_idx], i_pt, idx);
                        i_pt++;
                    }
                }
            } else{
                if((miss_dist[idx]<miss_dist[idx-1]) && (miss_dist[idx]<miss_dist[idx+1])){
                    local_min = miss_dist[idx];
                    if((i_pt < nPoints) && (local_min<4*r_c[idx]*r_c[idx])){
                        wake.rotor_wakes[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction[i_blade_idx].interaction_pts.wake_idx[i_pt] = idx;
                        wake.rotor_wakes[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction[i_blade_idx].interaction_pts.miss_dist[i_pt] = miss_dist[idx];
                        wake.rotor_wakes[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction[i_blade_idx].interaction_pts.r_c[i_pt] = r_c[idx];
                        fill_interactionPts(wake.rotor_wakes[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction[i_blade_idx], i_pt, idx);
                        i_pt++;
                    }
                }
            }
            
        }
    }   
} 

directionVec[] get_directionVec(string value, BWI)(auto ref BWI VortexInteraction) {
	
	immutable elements = VortexInteraction.BWI_inputs.length * chunk_size;

	directionVec[] d = new directionVec[elements];

	foreach(c_idx, ref chunk; VortexInteraction.BWI_inputs) {
		immutable out_start_idx = c_idx*chunk_size;

		immutable remaining = elements - out_start_idx;
		
		immutable out_end_idx = remaining > chunk_size ? (c_idx + 1)*chunk_size : out_start_idx + remaining;
		immutable in_end_idx = remaining > chunk_size ? chunk_size : remaining;

        foreach(idx; 0..in_end_idx){
            mixin("d[out_start_idx + idx].x = chunk."~value~"[idx].x;");
            mixin("d[out_start_idx + idx].y = chunk."~value~"[idx].y;");
            mixin("d[out_start_idx + idx].z = chunk."~value~"[idx].z;"); 
        }
		
	}
	return d;
}


void fill_interactionPts (BWI) (ref auto BWI interaction, size_t i_pt, size_t idx) {
    size_t[] bladeSec_idx;
    double[] gamma_w;
    directionVec[] r_blade;
    directionVec[] r_vortex;

    bladeSec_idx = get_blade_sec_idx!"bladeSec_idx"(interaction);
    gamma_w = get_BWIinputs!"gamma_w"(interaction);
    //r_blade = get_directionVec!"r_blade"(interaction);
    r_vortex = get_directionVec!"r_vortex"(interaction);

    interaction.interaction_pts.bladeSec_idx[i_pt] = bladeSec_idx[idx];
    interaction.interaction_pts.gamma_w[i_pt] = gamma_w[idx];
    //interaction.interaction_pts.r_blade[i_pt].x = r_blade[idx].x;
    //interaction.interaction_pts.r_blade[i_pt].y = r_blade[idx].y;
    //interaction.interaction_pts.r_blade[i_pt].z = r_blade[idx].z;
    interaction.interaction_pts.r_vortex[i_pt].x = r_vortex[idx].x;
    interaction.interaction_pts.r_vortex[i_pt].y = r_vortex[idx].y;
    interaction.interaction_pts.r_vortex[i_pt].z = r_vortex[idx].z;
}

void store_bladeSec_gamma (BWI, BS) (ref auto BWI tipVortex, BS blade_state){
    double[] gamma = get_wake_component!"gamma"(blade_state);
    double[] x = get_wake_component!"x"(blade_state);
    double[] y = get_wake_component!"y"(blade_state);
    double[] z = get_wake_component!"z"(blade_state);
    foreach(blade_idx; 0..tipVortex.length){
        foreach(i_pt; 0..nPoints){
            tipVortex[blade_idx].interaction_pts.gamma_sec[i_pt] = gamma[tipVortex[blade_idx].interaction_pts.bladeSec_idx[i_pt]];
            if(tipVortex[blade_idx].interaction_pts.bladeSec_idx[i_pt] != 0){
                tipVortex[blade_idx].interaction_pts.r_blade[i_pt].x = x[tipVortex[blade_idx].interaction_pts.bladeSec_idx[i_pt]] - x[tipVortex[blade_idx].interaction_pts.bladeSec_idx[i_pt]-1];
                tipVortex[blade_idx].interaction_pts.r_blade[i_pt].y = y[tipVortex[blade_idx].interaction_pts.bladeSec_idx[i_pt]] - y[tipVortex[blade_idx].interaction_pts.bladeSec_idx[i_pt]-1];
                tipVortex[blade_idx].interaction_pts.r_blade[i_pt].z = z[tipVortex[blade_idx].interaction_pts.bladeSec_idx[i_pt]] - z[tipVortex[blade_idx].interaction_pts.bladeSec_idx[i_pt]-1];
            }
        }
    }
}