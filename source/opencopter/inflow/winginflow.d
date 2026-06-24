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

import numd.linearalgebra.matrix;
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

alias WingInflow = WingInflowT!(ArrayContainer.none);

class WingInflowT(ArrayContainer AC = ArrayContainer.none) : InflowT!AC {

    //Frame* local_frame;

    alias WG = WingGeometryT!AC;
    alias WS = WingStateT!AC;
    alias WIS = WingInputStateT!AC;
    alias WLS = WingLiftSurfT!AC;

    private WG* wing;
    private WS* wing_state;
    private WIS* wing_input;
    private WLS* wing_lift_surf;

    private Chunk[] dC_L;
    private Chunk[] y;

    Frame* local_frame;
	Mat4 global_inverse;

    @nogc Frame* frame() {
        return local_frame;
    }

    @nogc Mat4 inverse_global_frame() {
		return global_inverse;
	}

    this(WG* _wing, WIS* _wing_input, WLS* _wing_lift_surf){

        _wing.frame.children ~= new Frame(Vec3(1, 0, 0), 0, Vec3(0, 0, 0), _wing.frame, _wing.frame.name ~ " inflow", "connection");
		local_frame = _wing.frame.children[$-1];

        wing = _wing;
        //wing_state = _wing_state;
        wing_input = _wing_input;
        wing_lift_surf = _wing_lift_surf;

        dC_L = new Chunk[2*wing.wing_parts[0].chunks.length];
        y = new Chunk[2*wing.wing_parts[0].chunks.length];

        foreach(ref _y; y) {
            _y[] = 0.0;
        }
        foreach(ref _dC_L; dC_L) {
            _dC_L[] = 0.0;
        }

        size_t c_idx = 0;
        foreach(p_idx, ref wing_part; wing.wing_parts) {
            //writeln("p_idx: ", p_idx);
            foreach(ref chunk; wing_part.ctrl_chunks) {
                //writeln("c_idx: ",c_idx);
                if(c_idx < wing_part.chunks.length){
                    y[c_idx] = chunk.ctrl_pt_y[];
                    c_idx++;
                }
                
            }
        }
    }
    
    void update(AircraftStateT!AC ac_state, WakeT!AC wake, double dt){
        bool calc_wake_ind_vel = true;
        bool wing_only_test = false;
        auto wing_state = ac_state.wing_states[].filter!(WS => WS.inflow_model.frame == this.frame).front;
        
        immutable num_span_chunks = wing_state.wing_part_states[0].chunks.length;
        immutable num_chord_pt =  wing_state.wing_part_states[0].ctrl_chunks.length/num_span_chunks;
        
        /*foreach(wp_idx, wp; wing_state.wing_part_states){
            foreach(ctrl_ch_idx, ctrl_chunk; wp.ctrl_chunks){
                writeln("wp_idx = ", wp_idx, "chunk_idx = ", ctrl_ch_idx, "\tu_t = ", ctrl_chunk.ctrl_pt_ut[]);
                writeln("wp_idx = ", wp_idx, "chunk_idx = ", ctrl_ch_idx, "\tu_p = ", ctrl_chunk.ctrl_pt_up[]);
                writeln("wp_idx = ", wp_idx, "chunk_idx = ", ctrl_ch_idx, "\taoa = ", ctrl_chunk.ctrl_pt_aoa[]);
            }
        }*/
        /*foreach(wp_idx, wp_surf; wing_lift_surf.wing_part_lift_surf){
            foreach(fl_idx, filament; wp_surf.spanwise_filaments){
                foreach(ch_idx, chunk; filament.chunks){
                    writeln("\nwp_idx = ", wp_idx, "\tfl_idx = ", fl_idx, "\tch_idx = ", ch_idx, "\tA_kl = ", chunk.A_kl[]);
                    writeln("\tgamma = ", chunk.gamma[]);
                }
            }
        }*/


        Chunk wing_aoa = wing_input.angle_of_attack; 
        debug writeln("going into for loop");
        foreach(wp_idx, wing_part; wing.wing_parts){
            auto combined_inflow = Vector!(4, Chunk)(0.0);
            foreach(ch_idx,ctrl_chunk; wing_part.ctrl_chunks){
                if (!wing_only_test){
                    auto ctrl_xyz = Vector!(4, Chunk)(1.0);
                    ctrl_xyz[0][] = ctrl_chunk.ctrl_pt_x[];
                    ctrl_xyz[1][] = ctrl_chunk.ctrl_pt_y[];
                    ctrl_xyz[2][] = ctrl_chunk.ctrl_pt_z[];

                    //writeln("\n wing_inflow");
                    //writeln("ctrl_chunk:", "\tx = ",  ctrl_chunk.ctrl_pt_x[]);
                    //writeln("ctrl_chunk:", "\ty = ",  ctrl_chunk.ctrl_pt_y[]);
                    //writeln("ctrl_chunk:", "\tz = ",  ctrl_chunk.ctrl_pt_z[]);

                    auto xyz_global = local_frame.global_matrix * ctrl_xyz; //back corrected!!

                    auto local_inflow = Vector!(4, Chunk)(0);

                    if(calc_wake_ind_vel){
                        //writeln("\n");
                    
                        immutable Chunk x = xyz_global[0][];
                        immutable Chunk y = xyz_global[1][];
                        immutable Chunk z = xyz_global[2][];

                        /*writeln("wing_global_x = ", x);
                        writeln("wing_global_y = ", y);
                        writeln("wing_global_z = ", z);*/

                        auto wake_vels = wake.compute_wake_induced_velocities_on_wing(ac_state, x, y, z);

                        local_inflow[0][] = wake_vels.v_x[];
                        local_inflow[1][] = wake_vels.v_y[];
                        local_inflow[2][] = wake_vels.v_z[];

                        /*writeln("wake induced velocities x on wing", local_inflow[0][]);
                        writeln("wake induced velocities y on wing", local_inflow[1][]);
                        writeln("wake induced velocities z on wing", local_inflow[2][]);*/

                        //writeln("\n");

                        combined_inflow += local_inflow;

                    }else{
                        foreach(rotor_idx, ref rotor; ac_state.rotor_states) {

                            auto xyz_tpp = rotor.inflow_model.frame.inverse_global_matrix * xyz_global;  //back_corrected!!         

                            immutable Chunk rotor_induced = rotor.inflow_model.inflow_at(xyz_tpp);   

                            auto local_inflow_due_to_rotor = Vector!(4, Chunk)(0);                     

					        local_inflow_due_to_rotor[2][] = rotor_induced[];

                            //writeln("local_inflow_due_to_rotor z =", local_inflow_due_to_rotor[2][]);

                            combined_inflow += rotor.inflow_model.frame.global_matrix * local_inflow_due_to_rotor; //back_corrected!!
                            //writeln("combined_inflow after rotor x =", combined_inflow[0][]);
                            //writeln("combined_inflow after rotor y =", combined_inflow[1][]);
                            //writeln("combined_inflow after rotor z =", combined_inflow[2][]);
                        }
                    }
                
                    //writeln("rotor inflows calculated", combined_inflow);
                    foreach(if_idx, ref wing; ac_state.wing_states){
                        if(wing.inflow_model != this) {
                            auto xyz_tpp = wing.inflow_model.frame.inverse_global_matrix * xyz_global; //back corrected!!

                            immutable Chunk wing_induced = wing.inflow_model.inflow_at(xyz_tpp);

                            auto local_inflow_due_to_wing = Vector!(4, Chunk)(0);

                            local_inflow_due_to_wing[2][] = wing_induced[];

                            combined_inflow += wing.inflow_model.frame.global_matrix * local_inflow_due_to_wing; //back_corrected!!
                            writeln("going into wing_inflow");
                        }
                    }
                }
                
                
                //writeln("wing inflows calculated", combined_inflow);

                combined_inflow[0][] += ac_state.freestream[0];
                combined_inflow[1][] += ac_state.freestream[1];
                combined_inflow[2][] += ac_state.freestream[2];                
                
                /*writeln("\n ctrl_chunk_idx = ", ch_idx, "ctrl_chunk = ", ctrl_chunk);
                writeln("combined_inflow_x = ", combined_inflow[0]);
                writeln("combined_inflow_y = ", combined_inflow[1]);
                writeln("combined_inflow_z = ", combined_inflow[2]);*/

                //auto inv_global_mat = local_frame.global_matrix.inverse.get();
                //writeln("inv_global_mat = ", inv_global_mat);
                immutable Vector!(4, Chunk) wing_local_inflow = local_frame.inverse_global_matrix*combined_inflow; //back corrected!!

                /*writeln("wing_local_inflow_x = ", wing_local_inflow[0]);
                writeln("wing_local_inflow_y = ", wing_local_inflow[1]);
                writeln("wing_local_inflow_z = ", wing_local_inflow[2]);*/


                wing_state.wing_part_states[wp_idx].ctrl_chunks[ch_idx].ctrl_pt_up[] = wing_local_inflow[2][];
                wing_state.wing_part_states[wp_idx].ctrl_chunks[ch_idx].ctrl_pt_ut[] = wing_local_inflow[0][];
                debug writeln("wp_idx = ", wp_idx, "\twing control point up = ",wing_state.wing_part_states[wp_idx].ctrl_chunks[ch_idx].ctrl_pt_up[]);
                debug writeln("wp_idx = ", wp_idx, "\twing control point ut = ",wing_state.wing_part_states[wp_idx].ctrl_chunks[ch_idx].ctrl_pt_ut[]);

                immutable wing_inflow_angle = atan2(wing_local_inflow[2], wing_local_inflow[0]);
                

                immutable Chunk effective_aoa = wing_inflow_angle[];
                debug writeln("effective_aoa = ", effective_aoa, "\twing_inflow_angle = ", wing_inflow_angle, "\n");
                
                wing_state.wing_part_states[wp_idx].ctrl_chunks[ch_idx].ctrl_pt_aoa[] = ctrl_chunk.camber[] - effective_aoa[];
                //writeln("camber = ", ctrl_chunk.camber, "effective aoa = ", effective_aoa);
                combined_inflow = Vector!(4, Chunk)(0.0);
            }
            //wing_state.wing_part_states[wp_idx].circulation_model.compute_wing_C_L(wing_wing_part, wing_state.wing_part_states[wp_idx]);
        }
        debug writeln("going into update_wing_circulation");
        update_wing_circulation(wing_state);
        update_wing_dC_L(wing_state);

        foreach(w_idx, ref ws; ac_state.wing_states) {
            if (ws.inflow_model.frame == this.frame) {
                compute_wing_C_L(wing_lift_surf,ws,wing);
            }
        }
	    
        //writeln("wing C_L (ac_state) =", ac_state.wing_states[0].C_L);   
        
    }

    void update_wing_circulation(WingStateT!AC wing_state){
        debug writeln("going into update_wing_circulation");

        debug writeln("wing_state.wing_part_states.length = ",wing_state.wing_part_states.length);
        debug writeln("wing_state.wing_part_states[0].chunks = ",wing_state.wing_part_states[0].chunks);
        debug writeln("wing_state.wing_part_states[0].chunks.length = ", wing_state.wing_part_states[0].chunks.length);
        debug writeln("wing_state.wing_part_states[0].ctrl_chunks.length/num_span_chunks = ", wing_state.wing_part_states[0].ctrl_chunks.length);
        immutable num_span_chunks = wing_state.wing_part_states[0].chunks.length;
        immutable num_chord_pt =  wing_state.wing_part_states[0].ctrl_chunks.length/num_span_chunks;

        
        foreach(wp_idx, wing_part; wing.wing_parts){
            foreach(chord_idx;0..num_chord_pt){
                foreach(ch_idx; 0..num_span_chunks){
                    /*immutable ctrl_ch_idx = chord_idx*num_span_chunks + ch_idx;
                    //writeln("ctrl_chunk_length = ",wing_state.wing_part_states[wp_idx].ctrl_chunks.length, "\tctrl_ch_idx = ", ctrl_ch_idx);
                    //writeln("wing control point ut = ",wing_state.wing_part_states[wp_idx].ctrl_chunks[ctrl_ch_idx].ctrl_pt_ut[]);
                    immutable Chunk u_t = wing_state.wing_part_states[wp_idx].ctrl_chunks[ctrl_ch_idx].ctrl_pt_ut;
                    immutable Chunk u_p = wing_state.wing_part_states[wp_idx].ctrl_chunks[ctrl_ch_idx].ctrl_pt_up;
                    immutable Chunk u_sqared = (u_t[]*u_t[] + u_p[]*u_p[]);
                    immutable Chunk u_inf = sqrt(u_sqared);*/
                    wing_state.wing_part_states[wp_idx].circulation_model.compute_d_gamma_coefficients(wing_lift_surf, wing_state.wing_part_states[wp_idx], wp_idx, ch_idx, chord_idx);
                }
            }

            foreach(chord_idx;0..num_chord_pt){
                foreach(ch_idx; 0..num_span_chunks){
                    immutable ctrl_ch_idx = chord_idx*num_span_chunks + ch_idx;
                    //writeln("ctrl_chunk_length = ",wing_state.wing_part_states[wp_idx].ctrl_chunks.length, "\tctrl_ch_idx = ", ctrl_ch_idx);
                    //writeln("wing control point ut = ",wing_state.wing_part_states[wp_idx].ctrl_chunks[ctrl_ch_idx].ctrl_pt_ut[]);
                    immutable Chunk u_t = wing_state.wing_part_states[wp_idx].ctrl_chunks[ctrl_ch_idx].ctrl_pt_ut;
                    immutable Chunk u_p = wing_state.wing_part_states[wp_idx].ctrl_chunks[ctrl_ch_idx].ctrl_pt_up;
                    immutable Chunk u_sqared = (u_t[]*u_t[] + u_p[]*u_p[]);
                    immutable Chunk u_inf = sqrt(u_sqared);
                    wing_state.wing_part_states[wp_idx].circulation_model.compute_bound_circulation(wing_lift_surf, wing_part, wp_idx, ch_idx, chord_idx, u_inf);
                }
            }
        }
    }

    void update_wing_dC_L(WingStateT!AC wing_state){
        immutable num_span_chunks = wing_state.wing_part_states[0].chunks.length;
        immutable num_chord_pt =  wing_state.wing_part_states[0].ctrl_chunks.length/num_span_chunks;
        double root_chord = wing.wing_parts[0].wing_root_chord; 
        size_t c_idx = 0;
        foreach(wp_idx, wing_part_state; wing_state.wing_part_states){
            foreach (span_idx; 0..num_span_chunks){
                wing_part_state.circulation_model.compute_dCl(wing_lift_surf,wing_part_state,wp_idx,span_idx);
                //wing_part_state.chunks[span_idx].dC_L[] /= root_chord; //circulation is multiplied by root chord in compute_d_gamma_circulation, so Cl need to be devided by it
                dC_L[c_idx][] = wing_part_state.chunks[span_idx].dC_L[];
                c_idx++;
            }

            //compute_wing_C_L(wing_lift_surf, wing, wing_state);
        }
    }

    /*void compute_wing_C_L(WLS, W, WS)(auto ref WLS wing_lift_surface, auto ref W wing, auto ref WS wing_state){

        //Calculating wing C_L using lifting line theory formula from the frink paper
        //Not 100% sure what how does the formula including l, k = 1 to N and M works for entire wing and not one wing part
        

        //wing_state.C_L = integrate_trapaziodal!"dC_L"(wing_state, wing);
        //debug writeln("wing C_L = ", wing_state.C_L);
    
        immutable AR = 2*(2*wing.wing_parts[0].wing_span)/(wing.wing_parts[0].wing_root_chord + wing.wing_parts[0].wing_tip_chord);
        immutable lamda = wing.wing_parts[0].wing_tip_chord/wing.wing_parts[0].wing_root_chord;
        immutable delta_m = 4*(lamda - 1)/(AR * (lamda +1));
        immutable M = wing.wing_parts[0].chunks.length* chunk_size; // number of spanwise nodes
        immutable N = wing_lift_surface.wing_part_lift_surf[0].spanwise_filaments.length ; //number chordwise vortex nodes

        //writeln("AR = ", AR);
        //writeln("lamda = ", lamda);
        //writeln("delta_m = ", delta_m);
        //writeln("M = ", M);
        //writeln("N = ", N);

        double term_1 = 0.0;
        double term_2 = 0.0;
        foreach(wp_idx, ref wp_lift_surface; wing_lift_surface.wing_part_lift_surf){
            foreach(fl_idx, ref filament; wp_lift_surface.spanwise_filaments){  //loop over all the filaments (chordwise)
                foreach(ch_idx, ref chunk; filament.chunks){ // loop over all the sapnwise chunks in a filament
                    foreach(c_idx; 0..chunk_size){
                        immutable k = fl_idx;
                        immutable l = ch_idx*chunk_size + c_idx;
                        
                        immutable theta_k = (2*k +1)*PI/(2*N);
                        immutable phi_l = (2*l +1)*PI/(2*M);
                        
                        term_1 += chunk.A_kl[c_idx]*sin(theta_k)*cos(phi_l);
                        term_2 += chunk.A_kl[c_idx]*sin(theta_k)*sin(phi_l)*sin(phi_l);

                        //writeln("wp_idx = ", wp_idx, "\tk = ", k, "\tl = ", l, "\ttheta_k = ", theta_k, "\tphi_l = ", phi_l, "\tA_kl = ", chunk.A_kl[c_idx], "\tterm_1 = ", term_1, "\tterm_2 = ", term_2);
                    }
                }                
            }
        }
        
        double C_L =  (AR*delta_m/8 + 1/(1+lamda)) * (wing_state.wing_part_states[0].chunks[0].dC_L[0] + PI*PI * term_1/(N*M)) + (delta_m*AR/16) * PI*PI*term_2/(N*M);
        debug writeln("C_L from formula = ", C_L);
        
        wing_state.C_L = C_L;            
    }*/

    void compute_wing_C_L(WLS, WS, W)(auto ref WLS wing_lift_surface, auto ref WS wing_state, auto ref W wing){

        double C_L = 0.0;
        immutable num_span_chunks = wing_state.wing_part_states[0].chunks.length;
        foreach(wp_idx, wing_part_state; wing_state.wing_part_states){
            foreach(span_idx; 0..wing_part_state.chunks.length){
                immutable Chunk dC_L = wing_part_state.chunks[span_idx].dC_L[];
                immutable Chunk span_y = wing.wing_parts[wp_idx].chunks[span_idx].y_span[];
                double y_1 = 0.0;
                double y_2 = 0.0;
                foreach(c1;0..chunk_size){
                    y_1 = wing_lift_surface.wing_part_lift_surf[wp_idx].spanwise_filaments[0].chunks[span_idx].y[c1];
                    if(c1 == 0){
                        if(span_idx == 0){
                        y_2 = 0;
                        }else{
                        y_2 = wing_lift_surface.wing_part_lift_surf[wp_idx].spanwise_filaments[0].chunks[span_idx-1].y[7];
                        }
                    }else{  
                        y_2 = wing_lift_surface.wing_part_lift_surf[wp_idx].spanwise_filaments[0].chunks[span_idx].y[c1 - 1];
                    }
                    double delta_y = y_1-y_2;
                    C_L += dC_L[c1]*delta_y/wing.wing_parts[wp_idx].wing_span;
                }
                
            }
        }

        wing_state.C_L = C_L;
    }

    InducedVelocities compute_wing_induced_vel_on_blade(immutable Chunk x, immutable Chunk y, immutable Chunk z){
        auto ind_vel = compute_wing_induced_vel(wing_lift_surf, x, y, z);

        return ind_vel;
    }
    
    Chunk inflow_at(immutable Chunk x, immutable Chunk y, immutable Chunk z, immutable Chunk x_e, double angle_of_attack){
        auto ind_vel = compute_wing_induced_vel(wing_lift_surf, x, y, z);

        immutable Chunk V_z = ind_vel.v_z[];
        return V_z;
    }

    Chunk inflow_at(immutable Vector!(4, Chunk) xyz) {
        auto ind_vel = compute_wing_induced_vel(wing_lift_surf, xyz[0], xyz[1], xyz[2]);
        immutable Chunk V_z = ind_vel.v_z[];
        return V_z;
    }

    double wake_skew() {
        return 0.0;
    }
}
