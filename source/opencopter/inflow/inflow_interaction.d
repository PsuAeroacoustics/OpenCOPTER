module inflow.inflow_interaction;

import opencopter.aircraft;
import opencopter.config;
import opencopter.inflow;
import opencopter.math;
import opencopter.math.blas;
import opencopter.math.lapacke;
import opencopter.memory;
import opencopter.wake;
import opencopter.vortexlattice;

import numd.calculus.integration.forwardeuler;
import numd.calculus.integration.rk4;
import numd.utility : linspace;

import std.array;
import std.algorithm;
import std.complex;
import std.conv;
import std.math : sin,cos,PI;
import std.range;
import std.stdio;
import std.typecons;

void get_ind_vel_on_rotor(RS,RG,RIS,WG,WIS,WS,I)(auto ref RS rotor_states,auto ref RG rotors, auto ref RIS rotor_inputs,auto ref WIS wing_inputs, auto ref WS wing_states, auto ref I[] inflows){
    
    size_t num_rotors = rotors.length; 
    double[num_rotors] v_z;
    double[num_rotors] v_x;
    Chunk x_rotor; Chunk y_rotor; Chunk z_rotor;
    Chunk psi_i = linspace(0,2*PI,8).to!double;

    foreach(rotor_idx, rotor; rotors){
        double cos_aoa = rotor_input[rotor_idx].cos_aoa;
        double sin_aoa = rotor_input[rotor_idx].sin_aoa;

        double half_rotor_radius = rotor.blades[0].chunks[$].r[$]/2;
        // check the relative distance of the x, y, zcoordinates for the rotors
        x_rotor[] = half_rotor_radius*cos(psi_i)[]*cos_aoa - rotor.origin[0];
        y_rotor[] = -(half_rotor_radius*sin(psi_i)[] + rotor.origin[1]);
        z_rotor[] = -half_rotor_radius*cos(psi_i)[]*sin_aoa + rotor.origin[2];
        foreach(inflow_idx,inflow; inflows){
            if(inflow_idx != rotor_idx){
                x_rotor[] = (x_rotor[] - rotor[inflow_idx].origin[0])*rotor_inputs[inflow_idx].cos_aoa;
                y_rotor[] = (-y_rotor[] - rotor[inflow_idx].origin[1]);
                z_rotor[] = (-z_rotor[] - rotor[inflow_idx].origin[2])*rotor_inputs[inflow_idx].sin_aoa;
                auto ind_vel = inflows.inflow_at(x_rotor,y_rotor,z_rotor,0,rotor_inpts[rotor_idx].angle_of_attack);
                if(inflow_idx < num_rotors){
                    v_z[rotor_idx] += -ind_vel.v_z[].mean*rotor_input[inflow_idx].cos_aoa;
                    v_x[rotor_idx] += ind_vel.v_z[].mean*rotor_input[inflow_idx].sin_aoa;
                } else {
                    v_z[rotor_idx] += ind_vel.v_z[].mean*wing_inputs[inflow_idx-num_rotors].cos_aoa;
                    v_x[rotor_idx] += -ind_vel.v_z[].mean*wing_inputs[inflow_idx-num_rotors].sin_aoa;
                }
                
                rotor_states[rotor_idx].advance_ratio = rotor_inputs[0].freestream_velocity + v_x[rotor_idx];
                rotor_states[rotor_idx].axial_advance_ratio = v_z[rotor_idx];  
            }
            
        }
    }

}