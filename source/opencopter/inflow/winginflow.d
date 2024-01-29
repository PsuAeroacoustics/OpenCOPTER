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

class WingInflowT(ArrayContainer AC = ArrayContainer.none) : InflowT!AC {

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

    void update(ref AircraftInputStateT!(AC) ac_input , ref AircraftT!(AC) aircraft, InflowT!AC[] inflows, double freestream_velocity, double advance_ratio, double axial_advance_ratio, double dt){
        update_impl(ac_input, aircraft, inflows, freestream_velocity, advance_ratio, axial_advance_ratio,dt);
    }

    void update(AircraftInputStateT!(AC)* ac_input , AircraftT!(AC)* aircraft, InflowT!AC[] inflows, double freestream_velocity, double advance_ratio, double axial_advance_ratio, double dt){
        update_impl(ac_input, aircraft, inflows, freestream_velocity, advance_ratio, axial_advance_ratio,dt);
    }

    void update_impl(AIS, AG)(auto ref AIS ac_input , auto ref AG aircraft, InflowT!AC[] inflows, double _freestream_velocity, double _advance_ratio, double _axial_advance_ratio, double dt){
        
        immutable num_span_chunks = wing_state.wing_part_states[0].chunks.length;
        immutable num_chord_pt =  wing_state.wing_part_states[0].ctrl_chunks.length/num_span_chunks;

        Chunk x; Chunk y; Chunk z;
        Chunk chunk_of_zeros = 0.0;
        Chunk rotor_inflow_x = chunk_of_zeros;
        Chunk rotor_inflow_z = chunk_of_zeros;
        Chunk x_wing_op; Chunk y_wing_op; Chunk z_wing_op;
        Chunk x_inflow; Chunk y_inflow; Chunk z_inflow;
        Vec3 freestream_velocity = [_freestream_velocity,0.0,0.0];

        Chunk wing_aoa = wing_input.angle_of_attack;     
        //writeln("wing_state_aoa = ", wing_input.angle_of_attack);   

        immutable num_rotors = ac_input.rotor_inputs.length;        
        foreach(wp_idx, wing_part; wing.wing_parts){
            //writeln("updating wing state properties");
            foreach(ch_idx,ctrl_chunk; wing_part.ctrl_chunks){
                //immutable span_idx = ch_idx%num_span_chunks;
                //immutable chord_idx = (ch_idx - ch_idx%num_span_chunks)/num_span_chunks;
                //x, y, z = conv_wing_to_OP_coord(ctrl_chunk.ctrl_pt_x,ctrl_chunk.ctrl_pt_y,ctrl_chunk.ctrl_pt_z);

                foreach(if_idx, inflow; inflows){
                    x_wing_op[] = -ctrl_chunk.ctrl_pt_x[]*wing_input.cos_aoa - ctrl_chunk.ctrl_pt_z[]*wing_input.sin_aoa + wing.origin[0];
                    y_wing_op[] = -ctrl_chunk.ctrl_pt_y[] + wing.origin[1];
                    z_wing_op[] = ctrl_chunk.ctrl_pt_z[]*wing_input.cos_aoa * ctrl_chunk.ctrl_pt_x[]*wing_input.sin_aoa + wing.origin[2];

                    //auto inflow_ind_vel = inflow.inflow_at(minus_ctrl_x[],minus_ctrl_y[],ctrl_chunk.ctrl_pt_z[],chunk_of_zeros,0);
                    //how to get rotor angle of attack for each rotor
                    if(if_idx < num_rotors){
                        x_wing_op[] = -ctrl_chunk.ctrl_pt_x[]*wing_input.cos_aoa - ctrl_chunk.ctrl_pt_z[]*wing_input.sin_aoa - aircraft.rotors[if_idx].origin[0];
                        y_wing_op[] = -ctrl_chunk.ctrl_pt_y[] - aircraft.rotors[if_idx].origin[1];
                        z_wing_op[] = ctrl_chunk.ctrl_pt_z[]*wing_input.cos_aoa * ctrl_chunk.ctrl_pt_x[]*wing_input.sin_aoa - aircraft.rotors[if_idx].origin[2];

                        x_inflow[] = x_wing_op[]*ac_input.rotor_inputs[if_idx].cos_aoa + z_wing_op[]*ac_input.rotor_inputs[if_idx].sin_aoa;
                        y_inflow[] = -y_wing_op[];
                        z_inflow[] = x_wing_op[]*ac_input.rotor_inputs[if_idx].sin_aoa + x_wing_op[]*ac_input.rotor_inputs[if_idx].cos_aoa;

                        //writeln("x_inflow_pt= ", x_inflow,"\ny_inflow_pt= ", y_inflow,"\nz_inflow_pt= ", z_inflow);

                        immutable Chunk inflow_ind_vel = inflow.inflow_at(x_inflow[], y_inflow[], z_inflow[], chunk_of_zeros, ac_input.rotor_inputs[if_idx].angle_of_attack);
                        //writeln("inflow_ind_vel_rotor: ", inflow_ind_vel);
                        //immutable cos_alpha = ac_input.rotor_inputs[if_idx].cos_aoa;
	                    //immutable sin_alpha = ac_input.rotor_inputs[if_idx].sin_aoa;
                        rotor_inflow_z[] += -inflow_ind_vel[] * ac_input.rotor_inputs[if_idx].cos_aoa * abs(ac_input.rotor_inputs[if_idx].angular_velocity * aircraft.rotors[if_idx].radius);
                        //rotor_inflow_y[] += -inflow_ind_vel.v_y[]*ac_input.rotor_inputs[if_idx].angular_velocity**ac_input.rotor_inputs.r_0[$];
                        rotor_inflow_x[] += -inflow_ind_vel[] * ac_input.rotor_inputs[if_idx].sin_aoa * abs(ac_input.rotor_inputs[if_idx].angular_velocity * aircraft.rotors[if_idx].radius);
                        //writeln(" rotor_inflow_x: ", rotor_inflow_x[], " rotor_inflow_z: ", rotor_inflow_z[]);
                    } else {
                        if( aircraft.wings[if_idx-num_rotors].origin != wing.origin){

                            //writeln("calculating influence due to wing");
                            x_wing_op[] = -ctrl_chunk.ctrl_pt_x[]*wing_input.cos_aoa - ctrl_chunk.ctrl_pt_z[]*wing_input.sin_aoa - aircraft.wings[if_idx - num_rotors].origin[0];
                            y_wing_op[] = -ctrl_chunk.ctrl_pt_y[] - aircraft.wings[if_idx - num_rotors].origin[1];
                            z_wing_op[] = ctrl_chunk.ctrl_pt_z[]*wing_input.cos_aoa * ctrl_chunk.ctrl_pt_x[]*wing_input.sin_aoa - aircraft.wings[if_idx - num_rotors].origin[2];

                            x_inflow[] = -x_wing_op[]*ac_input.wing_inputs[if_idx-num_rotors].cos_aoa - z_wing_op[]*ac_input.wing_inputs[if_idx - num_rotors].sin_aoa;
                            y_inflow[] = -y_wing_op[];
                            z_inflow[] = -x_wing_op[]*ac_input.wing_inputs[if_idx - num_rotors].sin_aoa + z_wing_op[]*ac_input.wing_inputs[if_idx - num_rotors].cos_aoa;
                            //immutable cos_alpha = cos(ac_input.wing_inputs[if_idx-num_rotors].aoa);
	                        //immutable sin_alpha = sin(ac_input.wing_inputs[if_idx-num_rotors].aoa);
                            immutable Chunk inflow_ind_vel = inflow.inflow_at(x_inflow[], y_inflow[], z_inflow[], chunk_of_zeros, ac_input.wing_inputs[if_idx - num_rotors].angle_of_attack);
                            //writeln("inflow_ind_vel_wing: ", inflow_ind_vel);
                            //writeln("inflow_ind_vel: ", inflow_ind_vel);
                            rotor_inflow_z[] += inflow_ind_vel[] * ac_input.wing_inputs[if_idx-num_rotors].cos_aoa;
                            //rotor_inflow_y[] += inflow_ind_vel.v_y[];
                            rotor_inflow_x[] += -inflow_ind_vel[] * ac_input.wing_inputs[if_idx-num_rotors].sin_aoa;
                        }                        
                    }
                }
                
                //writeln("rotor_inflow_x = ", rotor_inflow_x, "\trotor_inflow_z = ", rotor_inflow_z);
                immutable Chunk v_x = -(rotor_inflow_x[] - freestream_velocity[0]);
                immutable Chunk v_z = rotor_inflow_z[] + freestream_velocity[2];
                //writeln(" v_x: ", v_x[], " v_z: ", v_z[]);
                wing_state.wing_part_states[wp_idx].ctrl_chunks[ch_idx].ctrl_pt_up = rotor_inflow_z;
                wing_state.wing_part_states[wp_idx].ctrl_chunks[ch_idx].ctrl_pt_ut = v_x;
                immutable wing_inflow_angle = atan2(v_z,v_x);
                //writeln("v_x = ", v_x);
                //writeln("v_z = ", v_z);
                //writeln("freestream velocity = ", freestream_velocity[0]);
                //writeln("wing_aoa = ", wing_inflow_angle);
                immutable Chunk effective_aoa = wing_aoa[] + wing_inflow_angle[];
                //writeln("effective wing aoa = ", effective_aoa);
                //writeln("inflow angle: ", wing_inflow_angle[]);
                
                wing_state.wing_part_states[wp_idx].ctrl_chunks[ch_idx].ctrl_pt_aoa[] = ctrl_chunk.camber[] - effective_aoa[];
                //writeln("updating wing aoas");
                //writeln(" camber: ", ctrl_chunk.camber[], " effective aoa: ", effective_aoa[]);
                
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
                    //writeln("wing_ut = ",  wing_state.wing_part_states[wp_idx].ctrl_chunks[ctrl_ch_idx].ctrl_pt_ut);
                    wing_state.wing_part_states[wp_idx].circulation_model.compute_d_gamma_coefficients(wing_lift_surf, wing_state.wing_part_states[wp_idx], wp_idx, ch_idx, chord_idx, wing_state.wing_part_states[wp_idx].ctrl_chunks[ctrl_ch_idx].ctrl_pt_ut);
                    //writeln("wing_part = ", wp_idx, "\tchord_idx = ", chord_idx, "\tspan_chunk =", ch_idx ,"A_kl= ", wing_lift_surf.wing_part_lift_surf[wp_idx].spanwise_filaments[chord_idx].chunks[ch_idx].A_kl[]);
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

    void update_wing_dC_L(){
        immutable num_span_chunks = wing_state.wing_part_states[0].chunks.length;
        immutable num_chord_pt =  wing_state.wing_part_states[0].ctrl_chunks.length/num_span_chunks;
        double root_chord = wing.wing_parts[0].wing_root_chord; 
        foreach(wp_idx, wing_part_state; wing_state.wing_part_states){
            foreach (span_idx; 0..num_span_chunks){
                wing_part_state.circulation_model.compute_dCl(wing_lift_surf,wing_part_state,wp_idx,span_idx);
                wing_part_state.chunks[span_idx].dC_L[] /= root_chord; //circulation is multiplied by root chord in compute_d_gamma_circulation, so Cl need to be devided by it
            }
        }
    }

    InducedVelocities compute_wing_induced_vel_on_blade(immutable Chunk x, immutable Chunk y, immutable Chunk z){
        auto ind_vel = compute_wing_induced_vel(wing_lift_surf, x, y, z);

        return ind_vel;
    }    
    
    Chunk inflow_at(immutable Chunk x, immutable Chunk y, immutable Chunk z, immutable Chunk x_e, double angle_of_attack){
        auto ind_vel = compute_wing_induced_vel(wing_lift_surf, x, y, z);

        immutable Chunk neg_z = ind_vel.v_z[];
        //debug writeln("neg_z: ", neg_z);
        return neg_z;
    }
}