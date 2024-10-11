module opencopter.bwi;

import opencopter.aircraft;
import opencopter.config;
import opencopter.memory;
import opencopter.math;
import opencopter.wake;
import opencopter.bladeelement;
import opencopter.aircraft;

import std.math;
import std.container: DList ;
import std.range;

immutable size_t nPoints = 4;
immutable double BWI_factor = 4.0; // 2.0^2 as we are comparing distance^2

extern (C++) struct BWIinputsChunk {
    Chunk miss_dist;
    Chunk r_c_ave;
    Chunk gamma_w;
	Chunk gamma_sec;
    size_t[chunk_size] bladeSec_idx;
    Vec3[chunk_size] r_vortex;
}

extern (C++) struct InteractionPoints {
    size_t wake_idx;
    size_t bladeSec_idx;
    double miss_dist;
    double r_c;
    double gamma_w;
    double gamma_sec; 
    double[3] r_blade;
    double[3] r_vortex;
}

extern (C++) struct TipVortexInteractionT(ArrayContainer AC){
    import std.stdio : writeln;
    mixin ArrayDeclMixin!(AC, BWIinputsChunk, "BWI_inputs");

    auto interaction_pts = DList!InteractionPoints();
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

    size_t[] d;
    foreach (ref elem; VortexInteraction.interaction_pts) {
        mixin("d ~= elem."~value~";");
    }
    return d;
}

double[] get_interaction_point_components(string value, BWI)(auto ref BWI VortexInteraction){
    double[] d;
    foreach (ref elem; VortexInteraction.interaction_pts) {
        mixin("d ~= elem."~value~";");
    }
    return d;
}

double[][] get_interaction_point_directionVec(string value, BWI)(auto ref BWI VortexInteraction){
    Vec3[] vec;
    double[][] d;
    foreach (ref elem; VortexInteraction.interaction_pts) {
        mixin("d ~= elem."~value~";");
        
        /* foreach(idx; 0..2){
            mixin("d[][idx] ~= elem"~value~".[idx];");
        }  */
    }

    return d;
}


void calculate_BWI_points (W, BS)(auto ref W wake, auto ref BS blade_state, size_t rotor_idx, size_t blade_idx){
    // What we need
    // 1. Wake: x, y, z 
    // 2. r_c and miss distance
    // 3. local airfoil points: vehicle.aircraft.rotors.blades.chuncks (vtk.d ln340), it's calculated in bladeelement.d too!
    // Wake history is stored every 1 deg., therefore 1st interaction should be around 90 deg. 

    double[] miss_dist;
    double[] r_c;
    size_t i_pt=0;
    InteractionPoints interaction;

    foreach (i_blade_idx; 0..wake.rotor_wakes[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction.length){
        miss_dist = get_BWIinputs!"miss_dist"(wake.rotor_wakes[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction[i_blade_idx]);
        r_c = get_BWIinputs!"r_c_ave"(wake.rotor_wakes[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction[i_blade_idx]);
        auto bladeSec_idx = get_blade_sec_idx!"bladeSec_idx"(wake.rotor_wakes[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction[i_blade_idx]);
        auto gamma_w = get_BWIinputs!"gamma_w"(wake.rotor_wakes[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction[i_blade_idx]); 
        auto r_vortex = get_directionVec!"r_vortex"(wake.rotor_wakes[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction[i_blade_idx]);
        double[] gamma = get_wake_component!"gamma"(blade_state);
        double[] x = get_wake_component!"x"(blade_state);
        double[] y = get_wake_component!"y"(blade_state);
        double[] z = get_wake_component!"z"(blade_state);
        i_pt = 0;
        double local_min = miss_dist[0];
        foreach(idx; 0..miss_dist.length-1){
            if(idx == 0) {
                if(miss_dist[idx]<miss_dist[idx+1]){
                    local_min = miss_dist[idx];
                    if((i_pt < nPoints) && (local_min < BWI_factor*r_c[idx]*r_c[idx])){
                        interaction.wake_idx = idx;
                        interaction.miss_dist = miss_dist[idx];
                        interaction.r_c = r_c[idx];
                        interaction.bladeSec_idx = bladeSec_idx[idx];
                        interaction.gamma_w = gamma_w[idx];
                        interaction.r_vortex[0] = r_vortex[idx][0];
                        interaction.r_vortex[1] = r_vortex[idx][1];
                        interaction.r_vortex[2] = r_vortex[idx][2]; 
                        interaction.gamma_sec = gamma[bladeSec_idx[idx]];
                        interaction.r_blade[0] = x[bladeSec_idx[idx]]-x[bladeSec_idx[idx]-1];
                        interaction.r_blade[1] = y[bladeSec_idx[idx]]-y[bladeSec_idx[idx]-1];
                        interaction.r_blade[2] = z[bladeSec_idx[idx]]-z[bladeSec_idx[idx]-1];
                        i_pt++;
                        wake.rotor_wakes[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction[i_blade_idx].interaction_pts.insertBack(interaction);
                    }
                }
            } else{
                if((miss_dist[idx]<miss_dist[idx-1]) && (miss_dist[idx]<miss_dist[idx+1])){
                    local_min = miss_dist[idx];
                    if((i_pt < nPoints) && (local_min < BWI_factor*r_c[idx]*r_c[idx])){
                        interaction.wake_idx = idx;
                        interaction.miss_dist = miss_dist[idx];
                        interaction.r_c = r_c[idx];
                        interaction.bladeSec_idx = bladeSec_idx[idx];
                        interaction.gamma_w = gamma_w[idx];
                        interaction.r_vortex[0] = r_vortex[idx][0];
                        interaction.r_vortex[1] = r_vortex[idx][1];
                        interaction.r_vortex[2] = r_vortex[idx][2]; 
                        interaction.gamma_sec = gamma[bladeSec_idx[idx]];
                        interaction.r_blade[0] = x[bladeSec_idx[idx]]-x[bladeSec_idx[idx]-1];
                        interaction.r_blade[1] = y[bladeSec_idx[idx]]-y[bladeSec_idx[idx]-1];
                        interaction.r_blade[2] = z[bladeSec_idx[idx]]-z[bladeSec_idx[idx]-1]; 
                        i_pt++;
                        wake.rotor_wakes[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction[i_blade_idx].interaction_pts.insertBack(interaction);
                    }
                }
            }
            
        }
    }   
} 

Vec3[] get_directionVec(string value, BWI)(auto ref BWI VortexInteraction) {
	
	immutable elements = VortexInteraction.BWI_inputs.length * chunk_size;

	Vec3[] d = new Vec3[elements];

	foreach(c_idx, ref chunk; VortexInteraction.BWI_inputs) {
		immutable out_start_idx = c_idx*chunk_size;

		immutable remaining = elements - out_start_idx;
		
		//immutable out_end_idx = remaining > chunk_size ? (c_idx + 1)*chunk_size : out_start_idx + remaining;
		immutable in_end_idx = remaining > chunk_size ? chunk_size : remaining;

        foreach(idx; 0..in_end_idx){
            mixin("d[out_start_idx + idx] = chunk."~value~"[idx];");
        }
		
	} 
	return d;
} 
