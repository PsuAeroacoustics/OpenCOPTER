module opencopter.bwi;

import opencopter.aircraft;
import opencopter.config;
import opencopter.memory;
import opencopter.math;
import opencopter.wake;
import opencopter.bladeelement;
import opencopter.aircraft;

import std.math;
import std.container : DList;
import std.range;
import std.stdio : writeln;

import core.memory; // for GC.collect

immutable size_t nPoints = 6;
immutable double BWI_factor = 4.0; // 2.0^2 as we are comparing distance^2

extern (C++) struct BWIinputsChunk
{
    Chunk miss_dist;
    Chunk r_c_ave;
    Chunk gamma_w;
    Chunk gamma_sec;
    size_t[chunk_size] bladeSec_idx;
    Vec3[chunk_size] r_vortex;
    Vec3[chunk_size] pos_v;
    // Vec3[chunk_size] u_ind;
    Chunk dl;
}

extern (C++) struct InteractionPoints
{
    size_t wake_idx;
    size_t bladeSec_idx;
    double miss_dist;
    double r_c;
    double gamma_w;
    double gamma_sec;
    double[3] r_blade;
    double[3] r_vortex;
    double[3] r_blade_v;
    double secLen;
    double C_d;
    // double TKE;
    double l;
    double[3] normal;
}

extern (C++) struct TipVortexInteractionT(ArrayContainer AC)
{

    mixin ArrayDeclMixin!(AC, BWIinputsChunk, "BWI_inputs");

    DList!InteractionPoints interaction_pts;
    size_t length;

    this(size_t wake_history, size_t blade_elements)
    {

        assert(wake_history % chunk_size == 0);
        immutable num_chunks = wake_history / chunk_size;

        mixin(array_ctor_mixin!(AC, "BWIinputsChunk", "BWI_inputs", "num_chunks"));
        interaction_pts = DList!InteractionPoints();

        length = 0;
    }

}

extern (C++) struct VortexInteractionT(ArrayContainer AC)
{
    mixin ArrayDeclMixin!(AC, TipVortexInteractionT!(AC), "tip_vortex_interaction");

    size_t length;

    this(size_t num_blades, size_t wake_history, size_t blade_elements)
    {

        mixin(array_ctor_mixin!(AC, "TipVortexInteractionT!(AC)", "tip_vortex_interaction", "num_blades"));

        //writeln(mixin(array_ctor_mixin!(AC, "TipVortexInteractionT!(AC)", "tip_vortex_interaction", "num_blades")));

        foreach (ref interaction; tip_vortex_interaction)
        {
            interaction = TipVortexInteractionT!AC(wake_history, blade_elements);
        }

    }

}

struct VortexInteraction_multiRotorT(ArrayContainer AC)
{
    mixin ArrayDeclMixin!(AC, VortexInteractionT!(AC), "blade_vortex_interaction");

    this(size_t num_blades1, size_t num_blades2, size_t wake_history, size_t blade_elements)
    {
        //writeln( mixin(array_ctor_mixin!(AC, "VortexInteractionT!(AC)", "blade_vortex_interaction", "num_blades[idx]")));
        mixin(array_ctor_mixin!(AC, "VortexInteractionT!(AC)", "blade_vortex_interaction", "num_blades2"));

        //writeln(mixin(array_ctor_mixin!(AC, "TipVortexInteractionT!(AC)", "tip_vortex_interaction", "num_blades")));
        //writeln("blade_vortex_interaction.length:", blade_vortex_interaction.length);
        foreach (ref interaction; blade_vortex_interaction)
        {
            interaction = VortexInteractionT!AC(num_blades1, wake_history, blade_elements);
        }

    }

    /*this(size_t num_blades, size_t wake_history, size_t blade_elements) {
        
		mixin(array_ctor_mixin!(AC, "VortexInteractionT!(AC)", "blade_vortex_interaction", "num_blades"));

        //writeln(mixin(array_ctor_mixin!(AC, "TipVortexInteractionT!(AC)", "tip_vortex_interaction", "num_blades")));
        
        foreach(r_idx, ref interaction; blade_vortex_interaction) {
			interaction = VortexInteractionT!AC(num_blades, wake_history, blade_elements);
		}
        
	}*/

}

alias TipVortexInteraction = TipVortexInteractionT!(ArrayContainer.none);

double[] get_BWIinputs(string value, BWI)(auto ref BWI VortexInteraction)
{

    immutable elements = VortexInteraction.BWI_inputs.length * chunk_size;

    double[] d = new double[elements];

    foreach (c_idx, ref chunk; VortexInteraction.BWI_inputs)
    {
        immutable out_start_idx = c_idx * chunk_size;

        immutable remaining = elements - out_start_idx;

        immutable out_end_idx = remaining > chunk_size ? (c_idx + 1) * chunk_size
            : out_start_idx + remaining;
        immutable in_end_idx = remaining > chunk_size ? chunk_size : remaining;

        mixin("d[out_start_idx..out_end_idx] = chunk." ~ value ~ "[0..in_end_idx];");
    }
    return d;
}

void fill_BWIinputs(string value, BWI)(
    auto ref BWI VortexInteraction,
    ref double[] outBuf)
{
    immutable elements = VortexInteraction.BWI_inputs.length * chunk_size;
    if (outBuf.length < elements)
        outBuf.length = elements; // reuse allocation if possible

    foreach (c_idx, ref chunk; VortexInteraction.BWI_inputs)
    {
        immutable out_start_idx = c_idx * chunk_size;

        immutable remaining = elements - out_start_idx;

        immutable out_end_idx = remaining > chunk_size ? (c_idx + 1) * chunk_size
            : out_start_idx + remaining;
        immutable in_end_idx = remaining > chunk_size ? chunk_size : remaining;

        mixin("outBuf[out_start_idx..out_end_idx] = chunk." ~ value ~ "[0..in_end_idx];");
    }
}

size_t[] get_blade_sec_idx(string value, BWI)(auto ref BWI VortexInteraction)
{

    immutable elements = VortexInteraction.BWI_inputs.length * chunk_size;

    size_t[] d = new size_t[elements];

    foreach (c_idx, ref chunk; VortexInteraction.BWI_inputs)
    {
        immutable out_start_idx = c_idx * chunk_size;

        immutable remaining = elements - out_start_idx;

        immutable out_end_idx = remaining > chunk_size ? (c_idx + 1) * chunk_size
            : out_start_idx + remaining;
        immutable in_end_idx = remaining > chunk_size ? chunk_size : remaining;

        mixin("d[out_start_idx..out_end_idx] = chunk." ~ value ~ "[0..in_end_idx];");
    }
    return d;
}

size_t[] get_interaction_points(string value, BWI)(auto ref BWI VortexInteraction)
{

    size_t[] d;
    //assert(VortexInteraction !is null, "VortexInteraction is null");
    //assert(VortexInteraction.interaction_pts.length > 0, "interaction_pts is empty or uninitialized");
    //if (VortexInteraction is null || VortexInteraction.interaction_pts.length == 0){
    //  return [];
    //} else{
    foreach (ref elem; VortexInteraction.interaction_pts)
    {
        mixin("d ~= elem." ~ value ~ ";");
    }
    return d;
    //}
}

void fill_indexArray(string value, BWI)(
    auto ref BWI VortexInteraction,
    ref size_t[] outIdx)
{
    //immutable elements = VortexInteraction.interaction_pts.length();
    //if (outIdx.length < elements)
    //    outIdx.length = elements;

    foreach (ref elem; VortexInteraction.interaction_pts)
    {
        mixin("outIdx ~= elem." ~ value ~ ";");
    }
}

/*size_t[] get_interaction_points(string value, BWI)(auto ref BWI VortexInteraction)
{
    size_t[] d;

    //writeln("VortexInteraction.length:", VortexInteraction.length);

    // Check if interaction_pts is initialized and not empty
    static if (is(typeof(VortexInteraction.interaction_pts)))
    {
        if (VortexInteraction.interaction_pts.empty)
            return [];

        foreach (ref elem; VortexInteraction.interaction_pts)
        {
            mixin("d ~= elem." ~ value ~ ";");
        }
    }

    return d;
}*/

double[] get_interaction_point_components(string value, BWI)(auto ref BWI VortexInteraction)
{
    double[] d;
    foreach (ref elem; VortexInteraction.interaction_pts)
    {
        mixin("d ~= elem." ~ value ~ ";");
    }
    return d;
}

void fill_interaction_point_components(string value, BWI)(
    auto ref BWI VortexInteraction,
    ref double[] outArray)
{
    // immutable elements = VortexInteraction.interaction_pts.length;
    // if (outArray.length < elements)
    //     outArray.length = elements;

    foreach (ref elem; VortexInteraction.interaction_pts)
    {
        mixin("outArray ~= elem." ~ value ~ ";");
    }
}

double[][] get_interaction_point_directionVec(string value, BWI)(auto ref BWI VortexInteraction)
{
    Vec3[] vec;
    double[][] d;
    foreach (ref elem; VortexInteraction.interaction_pts)
    {
        mixin("d ~= elem." ~ value ~ ";");

        /* foreach(idx; 0..2){
            mixin("d[][idx] ~= elem"~value~".[idx];");
        }  */
    }

    return d;
}

void fill_interaction_point_directionVec(string value, BWI)(
    auto ref BWI VortexInteraction, double[][] outArray)
{
    // immutable elements = VortexInteraction.interaction_pts.length;
    // if (outArray.length < elements)
    //     outArray.length = elements;

    foreach (ref elem; VortexInteraction.interaction_pts)
    {
        mixin("outArray ~= elem." ~ value ~ ";");
    }
}

void calculate_BWI_points(W, BS, BG)(auto ref W wake, auto ref BS blade_state, size_t rotor_idx, size_t blade_idx, auto ref BG bladeGeom, double[3] normalVec, size_t iteration)
{
    // What we need
    // 1. Wake: x, y, z 
    // 2. r_c and miss distance
    // 3. local airfoil points: vehicle.aircraft.rotors.blades.chuncks (vtk.d ln340), it's calculated in bladeelement.d too!
    // Wake history is stored every 1 deg., therefore 1st interaction should be around 90 deg. 

    size_t i_pt = 0;
    InteractionPoints interaction;
    double[] C_d = get_wake_component!"dC_D"(blade_state);
    double[] r = get_wake_component!"r"(bladeGeom);
    double[] x = get_wake_component!"x"(blade_state);
    double[] y = get_wake_component!"y"(blade_state);
    double[] z = get_wake_component!"z"(blade_state);
    // debug writeln("x:", x);
    // debug writeln("y:", y);
    // debug writeln("z:", z);

    //GC.collect();
    //auto  stats1 = GC.stats();
    //writeln("1. GC used bytes: ", stats1.usedSize);
    // writeln("1. GC free bytes: ", stats1.freeSize);
    // writeln("1. GC total bytes: ", stats1.usedSize + stats1.freeSize);
    //writeln("1. GC page total: ", GC.Stats.pageSize);

    foreach (i_rotor_idx; 0 .. wake.rotor_wakes.length)
    {
        foreach (i_blade_idx; 0 .. wake.rotor_wakes[i_rotor_idx].interaction_perRotor[rotor_idx]
            .blade_vortex_interaction[blade_idx].tip_vortex_interaction.length)
        {
            wake.rotor_wakes[i_rotor_idx].interaction_perRotor[rotor_idx].blade_vortex_interaction[blade_idx]
                .tip_vortex_interaction[i_blade_idx].interaction_pts.clear();

            // GC.collect();
            // auto  stats2 = GC.stats();
            // writeln("iteration:", iteration, "2. GC used bytes: ", stats2.usedSize);
            // writeln("rotor_idx:", rotor_idx, "blade_idx:", blade_idx, "i_rotor_idx:", i_rotor_idx," i_blade_idx:", i_blade_idx);
            // writeln("2. GC free bytes: ", stats2.freeSize);
            // writeln("2. GC total bytes: ", stats2.usedSize + stats2.freeSize);
            //writeln("2. GC page total: ", GC.Stats.pageSize);

            double[] miss_dist;
            double[] r_c;

            miss_dist = get_BWIinputs!"miss_dist"(
                wake.rotor_wakes[i_rotor_idx].interaction_perRotor[rotor_idx].blade_vortex_interaction[blade_idx]
                    .tip_vortex_interaction[i_blade_idx]);
            r_c = get_BWIinputs!"r_c_ave"(
                wake.rotor_wakes[i_rotor_idx].interaction_perRotor[rotor_idx].blade_vortex_interaction[blade_idx]
                    .tip_vortex_interaction[i_blade_idx]);
            auto bladeSec_idx = get_blade_sec_idx!"bladeSec_idx"(
                wake.rotor_wakes[i_rotor_idx].interaction_perRotor[rotor_idx].blade_vortex_interaction[blade_idx]
                    .tip_vortex_interaction[i_blade_idx]);
            //writeln("bladeSec_idx:", bladeSec_idx);
            auto gamma_w = get_BWIinputs!"gamma_w"(
                wake.rotor_wakes[i_rotor_idx].interaction_perRotor[rotor_idx].blade_vortex_interaction[blade_idx]
                    .tip_vortex_interaction[i_blade_idx]);
            auto r_vortex = get_directionVec!"r_vortex"(
                wake.rotor_wakes[i_rotor_idx].interaction_perRotor[rotor_idx].blade_vortex_interaction[blade_idx]
                    .tip_vortex_interaction[i_blade_idx]);
            auto r_b = get_directionVec!"pos_v"(
                wake.rotor_wakes[i_rotor_idx].interaction_perRotor[rotor_idx].blade_vortex_interaction[blade_idx]
                    .tip_vortex_interaction[i_blade_idx]);
            double[] gamma = get_wake_component!"gamma"(blade_state);

            double[] dl = get_BWIinputs!"dl"(
                wake.rotor_wakes[i_rotor_idx].interaction_perRotor[rotor_idx].blade_vortex_interaction[blade_idx]
                    .tip_vortex_interaction[i_blade_idx]);

            // double[] dl;
            // fill_BWIinputs!"dl"(
            //     wake.rotor_wakes[i_rotor_idx].interaction_perRotor[rotor_idx]
            //         .blade_vortex_interaction[blade_idx]
            //         .tip_vortex_interaction[i_blade_idx], dl);

            double l = 0.0;
            //auto u_ind = get_directionVec!"u_ind"(wake.rotor_wakes[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction[i_blade_idx]);
            // double[] miss_dist;
            // double[] r_c;
            // fill_BWIinputs!"miss_dist"(
            //     wake.rotor_wakes[i_rotor_idx].interaction_perRotor[rotor_idx]
            //         .blade_vortex_interaction[blade_idx]
            //         .tip_vortex_interaction[i_blade_idx], miss_dist);

            // fill_BWIinputs!"r_c_ave"(
            //     wake.rotor_wakes[i_rotor_idx].interaction_perRotor[rotor_idx]
            //         .blade_vortex_interaction[blade_idx]
            //         .tip_vortex_interaction[i_blade_idx], r_c);
            // auto bladeSec_idx = get_blade_sec_idx!"bladeSec_idx"(
            //     wake.rotor_wakes[i_rotor_idx].interaction_perRotor[rotor_idx]
            //         .blade_vortex_interaction[blade_idx]
            //         .tip_vortex_interaction[i_blade_idx]);
            // //writeln("bladeSec_idx:", bladeSec_idx);
            // double[] gamma_w;
            // fill_BWIinputs!"gamma_w"(
            //     wake.rotor_wakes[i_rotor_idx].interaction_perRotor[rotor_idx]
            //         .blade_vortex_interaction[blade_idx]
            //         .tip_vortex_interaction[i_blade_idx], gamma_w);
            // Vec3[] r_vortex;

            // fill_directionVec!"r_vortex"(
            //     wake.rotor_wakes[i_rotor_idx].interaction_perRotor[rotor_idx]
            //         .blade_vortex_interaction[blade_idx]
            //         .tip_vortex_interaction[i_blade_idx], r_vortex);
            // Vec3[] r_b;
            // fill_directionVec!"pos_v"(
            //     wake.rotor_wakes[i_rotor_idx].interaction_perRotor[rotor_idx]
            //         .blade_vortex_interaction[blade_idx]
            //         .tip_vortex_interaction[i_blade_idx], r_b);
            // double[] gamma = get_wake_component!"gamma"(blade_state);

            i_pt = 0;
            double local_min = miss_dist[0];
            foreach (idx; 0 .. miss_dist.length - 1)
            {
                if (idx == 0)
                {
                    l = 0.0;
                    if (miss_dist[idx] < miss_dist[idx + 1])
                    {
                        local_min = miss_dist[idx];
                        if ((i_pt < nPoints) && (local_min < BWI_factor * r_c[idx] * r_c[idx]))
                        {
                            interaction.wake_idx = idx;
                            interaction.miss_dist = miss_dist[idx];
                            interaction.r_c = r_c[idx];
                            interaction.bladeSec_idx = bladeSec_idx[idx];
                            interaction.gamma_w = gamma_w[idx];
                            interaction.r_vortex[0] = r_vortex[idx][0];
                            interaction.r_vortex[1] = r_vortex[idx][1];
                            interaction.r_vortex[2] = r_vortex[idx][2];
                            interaction.gamma_sec = gamma[bladeSec_idx[idx]];
                            if (bladeSec_idx[idx] == 0)
                            {
                                interaction.r_blade[0] = x[bladeSec_idx[idx]] - x[bladeSec_idx[idx] + 1];
                                interaction.r_blade[1] = y[bladeSec_idx[idx]] - y[bladeSec_idx[idx] + 1];
                                interaction.r_blade[2] = z[bladeSec_idx[idx]] - z[bladeSec_idx[idx] + 1];
                            }
                            else
                            {
                                interaction.r_blade[0] = x[bladeSec_idx[idx]] - x[bladeSec_idx[idx] - 1];
                                interaction.r_blade[1] = y[bladeSec_idx[idx]] - y[bladeSec_idx[idx] - 1];
                                interaction.r_blade[2] = z[bladeSec_idx[idx]] - z[bladeSec_idx[idx] - 1];
                            }
                            interaction.r_blade_v[0] = r_b[idx][0] - x[bladeSec_idx[idx]];
                            interaction.r_blade_v[1] = r_b[idx][1] - y[bladeSec_idx[idx]];
                            interaction.r_blade_v[2] = r_b[idx][2] - z[bladeSec_idx[idx]];
                            interaction.secLen = (x[bladeSec_idx[idx]] - x[0]) * (
                                x[bladeSec_idx[idx]] - x[0]) + (y[bladeSec_idx[idx]] - y[0]) * (
                                y[bladeSec_idx[idx]] - y[0]) + (z[bladeSec_idx[idx]] - z[0]) * (
                                z[bladeSec_idx[idx]] - z[0]);
                            interaction.l = l;
                            if (bladeSec_idx[idx] == 0)
                            {
                                interaction.C_d = C_d[0] / (r[1] - r[0]);
                            }
                            else
                            {
                                double dr = r[bladeSec_idx[idx]] - r[bladeSec_idx[idx] - 1];
                                interaction.C_d = 0.5 * (
                                    C_d[bladeSec_idx[idx]] + C_d[bladeSec_idx[idx] - 1]) / dr;
                            }
                            //interaction.C_d = sumC_d(C_d); 
                            interaction.normal = normalVec;
                            i_pt++;
                            wake.rotor_wakes[i_rotor_idx].interaction_perRotor[rotor_idx].blade_vortex_interaction[blade_idx]
                                .tip_vortex_interaction[i_blade_idx].interaction_pts.insertBack(
                                    interaction);
                        }
                    }
                }
                else
                {
                    l = l + dl[idx - 1];
                    if ((miss_dist[idx] < miss_dist[idx - 1]) && (
                            miss_dist[idx] < miss_dist[idx + 1]))
                    {
                        local_min = miss_dist[idx];
                        if ((i_pt < nPoints) && (local_min < BWI_factor * r_c[idx] * r_c[idx]))
                        {
                            interaction.wake_idx = idx;
                            interaction.miss_dist = miss_dist[idx];
                            interaction.r_c = r_c[idx];
                            interaction.bladeSec_idx = bladeSec_idx[idx];
                            //writeln("bladeSec_idx[idx]:", bladeSec_idx[idx]);
                            interaction.gamma_w = gamma_w[idx];
                            interaction.r_vortex[0] = r_vortex[idx][0];
                            interaction.r_vortex[1] = r_vortex[idx][1];
                            interaction.r_vortex[2] = r_vortex[idx][2];
                            interaction.gamma_sec = gamma[bladeSec_idx[idx]];
                            if (bladeSec_idx[idx] == 0)
                            {
                                interaction.r_blade[0] = x[bladeSec_idx[idx]] - x[bladeSec_idx[idx] + 1];
                                interaction.r_blade[1] = y[bladeSec_idx[idx]] - y[bladeSec_idx[idx] + 1];
                                interaction.r_blade[2] = z[bladeSec_idx[idx]] - z[bladeSec_idx[idx] + 1];
                            }
                            else
                            {
                                interaction.r_blade[0] = x[bladeSec_idx[idx]] - x[bladeSec_idx[idx] - 1];
                                interaction.r_blade[1] = y[bladeSec_idx[idx]] - y[bladeSec_idx[idx] - 1];
                                interaction.r_blade[2] = z[bladeSec_idx[idx]] - z[bladeSec_idx[idx] - 1];
                            }
                            interaction.r_blade_v[0] = r_b[idx][0] - x[bladeSec_idx[idx]];
                            interaction.r_blade_v[1] = r_b[idx][1] - y[bladeSec_idx[idx]];
                            interaction.r_blade_v[2] = r_b[idx][2] - z[bladeSec_idx[idx]];
                            interaction.secLen = (x[bladeSec_idx[idx]] - x[0]) * (
                                x[bladeSec_idx[idx]] - x[0]) + (y[bladeSec_idx[idx]] - y[0]) * (
                                y[bladeSec_idx[idx]] - y[0]) + (z[bladeSec_idx[idx]] - z[0]) * (
                                z[bladeSec_idx[idx]] - z[0]);
                            interaction.l = l;
                            if (bladeSec_idx[idx] == 0)
                            {
                                interaction.C_d = C_d[0] / (r[1] - r[0]);
                            }
                            else
                            {
                                double dr = r[bladeSec_idx[idx]] - r[bladeSec_idx[idx] - 1];
                                interaction.C_d = 0.5 * (
                                    C_d[bladeSec_idx[idx]] + C_d[bladeSec_idx[idx] - 1]) / dr;
                            }
                            //interaction.C_d = sumC_d(C_d);  
                            interaction.normal = normalVec;
                            i_pt++;
                            wake.rotor_wakes[i_rotor_idx].interaction_perRotor[rotor_idx].blade_vortex_interaction[blade_idx]
                                .tip_vortex_interaction[i_blade_idx].interaction_pts.insertBack(
                                    interaction);
                            //writeln("i_rotor_idx:", i_rotor_idx, "rotor_idx:", rotor_idx, "blade_idx", blade_idx, "i_blade_idx", i_blade_idx);
                            //writeln("2. bwi, interaction point added ");
                        }
                    }
                }

            }
        }
    }
    //GC.collect();
    //auto  stats3 = GC.stats();
    //writeln(rotor_idx, blade_idx, "iteration:", iteration, "3. GC used bytes: ", stats3.usedSize);
    // writeln("3. GC free bytes: ", stats3.freeSize);
    // writeln("3. GC total bytes: ", stats3.usedSize + stats3.freeSize);
    //writeln("3. GC page total: ", GC.Stats.pageSize);
}

Vec3[] get_directionVec(string value, BWI)(auto ref BWI VortexInteraction)
{

    immutable elements = VortexInteraction.BWI_inputs.length * chunk_size;

    Vec3[] d = new Vec3[elements];

    foreach (c_idx, ref chunk; VortexInteraction.BWI_inputs)
    {
        immutable out_start_idx = c_idx * chunk_size;

        immutable remaining = elements - out_start_idx;

        //immutable out_end_idx = remaining > chunk_size ? (c_idx + 1)*chunk_size : out_start_idx + remaining;
        immutable in_end_idx = remaining > chunk_size ? chunk_size : remaining;

        foreach (idx; 0 .. in_end_idx)
        {
            mixin("d[out_start_idx + idx] = chunk." ~ value ~ "[idx];");
        }

    }
    return d;
}

void fill_directionVec(string value, BWI)(
    auto ref BWI VortexInteraction,
    ref Vec3[] outVec)
{
    immutable elements = VortexInteraction.BWI_inputs.length;
    if (outVec.length < elements)
        outVec.length = elements;

    foreach (c_idx, ref chunk; VortexInteraction.BWI_inputs)
    {
        immutable out_start_idx = c_idx * chunk_size;

        immutable remaining = elements - out_start_idx;

        //immutable out_end_idx = remaining > chunk_size ? (c_idx + 1)*chunk_size : out_start_idx + remaining;
        immutable in_end_idx = remaining > chunk_size ? chunk_size : remaining;

        foreach (idx; 0 .. in_end_idx)
        {
            mixin("outVec[out_start_idx + idx] = chunk." ~ value ~ "[idx];");
        }

    }
}

double sumC_d(double[] arr)
{
    double result = 0.0;
    foreach (elem; arr)
    {
        result += elem;
    }
    return result;
}
