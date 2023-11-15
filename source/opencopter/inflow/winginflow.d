module opencopter.inflow.winginflow;

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
import std.math;
import std.range;
import std.stdio;
import std.typecons;

/*Chunk conv_wing_to_OP_coord(Chunk x_data, Chunk y_data, Chunk z_data){
    immutable Chunk multiply_by_minus_one = -1.0;
    Chunk x_data_out = multiply_by_minus_one[]*x_data[];
    Chunk y_data_out = multiply_by_minus_one[]*y_data[];
    Chunk z_data_out = z_data; 
    //immutable Chunk z_data_tmp = multiply_by_minus_one[]*z_data;
    return x_data_out,y_data_out,z_data_out;
}*/

alias WingInflow = WingInflowT!(ArrayContainer.none);

class WingInflowT(ArrayContainer AC = ArrayContainer.none) : Inflow {

    alias WG = WingGeometryT!AC;
    alias WS = WingStateT!AC;
    alias WIS = WingInputStateT!AC;
    alias WLS = WingLiftSurfT!AC;

    private WG* wing;
    private WS* wing_state;
    private WIS* wing_input;
    private WLS* wing_lift_surf;

    this(WG* _wing, WS* _wing_state, WIS* _wing_input, WLS* _wing_lift_surf){
        wing= _wing;
        wing_state = _wing_state;
        wing_input = _wing_input;
        wing_lift_surf = _wing_lift_surf;
    }

    private void update(ArrayContainer AC = ArrayContainer.None)(ref AircraftInputStateT!AC ac_input , Inflow[] inflows, double freestream_velocity, double advance_ratio, double axial_advance_ratio, double dt){
        update_impl(ac_input, inflows, freestream_velocity, advance_ratio, axial_advance_ratio,dt);
    }

    private void update(ArrayContainer AC = ArrayContainer.None)(AircraftInputStateT!AC* ac_input , Inflow[] inflows, double freestream_velocity, double advance_ratio, double axial_advance_ratio, double dt){
        update_impl(*ac_input, inflows, freestream_velocity, advance_ratio, axial_advance_ratio,dt);
    }

    void update_impl(ArrayContainer AC = ArrayContainer.None)(ref AircraftInputStateT!AC ac_input , Inflow[] inflows, double _freestream_velocity, double _advance_ratio, double _axial_advance_ratio, double dt){
        
        immutable num_span_chunks = wing_state.wing_part_states[0].chunks.length;
        immutable num_chord_pt =  wing_state.wing_part_states[0].ctrl_chunks.length/num_span_chunks;

        Chunk x; Chunk y; Chunk z;
        Chunk rotor_inflow_x; Chunk rotor_inflow_y; Chunk rotor_inflow_z;

        Vec3 freestream_velocity = [_freestream_velocity,0.0,0.0];

        Chunk wing_aoa = wing_input.angle_of_attack;        

        immutable num_rotors = ac_input.rotor_inputs.length;        
        foreach(wp_idx, wing_part; wing.wing_parts){
            foreach(ch_idx,ctrl_chunk; wing_part.ctrl_chunks){
                //immutable span_idx = ch_idx%num_span_chunks;
                //immutable chord_idx = (ch_idx - ch_idx%num_span_chunks)/num_span_chunks;
                //x, y, z = conv_wing_to_OP_coord(ctrl_chunk.ctrl_pt_x,ctrl_chunk.ctrl_pt_y,ctrl_chunk.ctrl_pt_z);
                foreach(if_idx, inflow; inflows){
                    auto inflow_ind_vel = inflow.inflow_at(-ctrl_chunk.ctrl_pt_x,-ctrl_chunk.ctrl_pt_y,ctrl_chunk.ctrl_pt_z,0,0);
                    //how to get rotor angle of attack for each rotor
                    if(if_idx >= num_rotors){
                        immutable cos_alpha = ac_input.rotor_inputs[if_idx].cos_aoa;
	                    immutable sin_alpha = ac_input.rotor_inputs[if_idx].sin_aoa;
                        rotor_inflow_z[] += (-inflow_ind_vel.v_x[]*sin_alpha - inflow_ind_vel.v_z[]*cos_alpha)*ac_input.rotor_inputs[if_idx].angular_velocity*ac_input.rotor_inputs.r_0[$];
                        rotor_inflow_y[] += -inflow_ind_vel.v_y[]*ac_input.rotor_inputs[if_idx].angular_velocity**ac_input.rotor_inputs.r_0[$];
                        rotor_inflow_x[] += (inflow_ind_vel.v_x[]*cos_alpha - inflow_ind_vel.v_z[].sin_alpha)*ac_input.rotor_inputs[if_idx].angular_velocity*ac_input.rotor_inputs.r_0[$];
                    } else {
                        immutable cos_alpha = cos(ac_input.wing_inputs[if_idx-num_rotors].aoa);
	                    immutable sin_alpha = sin(ac_input.wing_inputs[if_idx-num_rotors].aoa);
                        rotor_inflow_z[] += inflow_ind_vel.v_x[]*sin_alpha + inflow_ind_vel.v_z[]*cos_alpha;
                        rotor_inflow_y[] += inflow_ind_vel.v_y[];
                        rotor_inflow_x[] += inflow_ind_vel.v_x[]*cos_alpha - inflow_ind_vel.v_z[].sin_alpha;
                    }
                }
                
                immutable Chunk v_x = -rotor_inflow_x[] + freestream_velocity[0];
                immutable Chunk v_z = -rotor_inflow_z[] + freestream_velocity[2];
                wing_state.wing_part_states[wp_idx].ctrl_chunks[ch_idx].ctrl_pt_up = rotor_inflow_z;
                wing_state.wing_part_states[wp_idx].ctrl_chunks[ch_idx].ctrl_pt_ut = v_x;
                immutable wing_inflow_angle = atan2(v_z,v_x);
                immutable Chunk effective_aoa = wing_aoa[] - wing_inflow_angle[];
                
                wing_state.wing_part_states[wp_idx].ctrl_chunks[ch_idx].ctrl_pt_aoa = ctrl_chunk.camber[] - effective_aoa[];
            }
        }
    }

    void update_wing_circulation(){
        immutable num_span_chunks = wing_state.wing_part_states[0].chunks.length;
        immutable num_chord_pt =  wing_state.wing_part_states[0].ctrl_chunks.length/num_span_chunks;

        foreach(wp_idx, wing_part; wing.wing_parts){
            foreach(chord_idx;0..num_chord_pt){
                foreach(ch_idx; 0..num_span_chunks){
                    immutable ctrl_ch_idx = chord_idx*num_span_chunks + ch_idx;
                    wing_state.wing_part_states[wp_idx].circulation_model.compute_d_gamma_coefficients(wing_lift_surf, wing_state.wing_part_states[wp_idx], wp_idx, ch_idx, chord_idx, wing_state.wing_part_states[wp_idx].ctrl_chunks[ctrl_ch_idx].ctrl_pt_ut);
                }
            }

            foreach(chord_idx;0..num_chord_pt){
                foreach(ch_idx; 0..num_span_chunks){
                    //immutable ctrl_ch_idx = chord_idx*num_span_chunks + ch_idx;                    
                    wing_state.wing_part_states[wp_idx].circulation_model.compute_bound_circulation(wing_lift_surf, wing_part, wp_idx, ch_idx, chord_idx);
                } 
            }
        }
    }
    
    
    Chunk inflow_at(immutable Chunk x, immutable Chunk y, immutable Chunk z, immutable Chunk x_e, double angle_of_attack){
        auto ind_vel = compute_wing_induced_vel(wing_lift_surf, x, y, z);

        immutable Chunk neg_z = -ind_vel.v_z[];
        debug writeln("neg_z: ", neg_z);
        return neg_z;
    }
}