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

    Frame* inflow_frame;

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

    @nogc Frame* frame() {
        return inflow_frame;
    }

    this(WG* _wing, WS* _wing_state, WIS* _wing_input, WLS* _wing_lift_surf){
        wing = _wing;
        wing_state = _wing_state;
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
            foreach(ref chunk; wing_part.ctrl_chunks) {
                y[c_idx] = chunk.ctrl_pt_y[];
                c_idx++;
            }
        }
    }
    
    void update(AircraftStateT!AC ac_state, double dt) {
        
        immutable num_span_chunks = wing_state.wing_part_states[0].chunks.length;
        immutable num_chord_pt =  wing_state.wing_part_states[0].ctrl_chunks.length/num_span_chunks;

        auto combined_inflow = Vector!(4, Chunk)(0.0);

        Chunk wing_aoa = wing_input.angle_of_attack;     

        foreach(wp_idx, wing_part; wing.wing_parts){
            foreach(ch_idx,ctrl_chunk; wing_part.ctrl_chunks) {
                auto ctrl_xyz = Vector!(4, Chunk)(1.0);
                ctrl_xyz[0][] = ctrl_chunk.ctrl_pt_x[];
                ctrl_xyz[1][] = ctrl_chunk.ctrl_pt_y[];
                ctrl_xyz[2][] = ctrl_chunk.ctrl_pt_z[];

                foreach(if_idx, ref rotor; ac_state.rotor_states) {
                    auto xyz_tpp = inflow_frame.inverse_global_matrix * ctrl_xyz;

                    immutable Chunk rotor_induced = rotor.inflow_model.inflow_at(xyz_tpp);

                    auto local_inflow = Vector!(4, Chunk)(0);

					local_inflow[2][] = rotor_induced[];

                    combined_inflow += rotor.inflow_model.frame.global_matrix * local_inflow;
                }

                foreach(if_idx, ref wing; ac_state.wing_states) {
                    if(wing.inflow_model != this) {
                        auto xyz_tpp = inflow_frame.inverse_global_matrix * ctrl_xyz;

                        immutable Chunk wing_induced = wing.inflow_model.inflow_at(xyz_tpp);

                        auto local_inflow = Vector!(4, Chunk)(0);

                        local_inflow[2][] = wing_induced[];

                        combined_inflow += wing.inflow_model.frame.global_matrix * local_inflow;
                    }
                }
                
                combined_inflow[0][] += ac_state.freestream[0];
                combined_inflow[1][] += ac_state.freestream[1];
                combined_inflow[2][] += ac_state.freestream[2];

                immutable Vector!(4, Chunk) wing_local_inflow = inflow_frame.global_matrix*combined_inflow;

                wing_state.wing_part_states[wp_idx].ctrl_chunks[ch_idx].ctrl_pt_up[] = wing_local_inflow[2][];
                wing_state.wing_part_states[wp_idx].ctrl_chunks[ch_idx].ctrl_pt_ut[] = wing_local_inflow[0][];

                immutable wing_inflow_angle = atan2(wing_local_inflow[2], wing_local_inflow[0]);

                immutable Chunk effective_aoa = wing_aoa[] + wing_inflow_angle[];
                
                wing_state.wing_part_states[wp_idx].ctrl_chunks[ch_idx].ctrl_pt_aoa[] = ctrl_chunk.camber[] - effective_aoa[];
            }
        }

        update_wing_circulation();
	    update_wing_dC_L();
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
                    wing_state.wing_part_states[wp_idx].circulation_model.compute_bound_circulation(wing_lift_surf, wing_part, wp_idx, ch_idx, chord_idx);
                }
            }
        }
    }

    void update_wing_dC_L(){
        immutable num_span_chunks = wing_state.wing_part_states[0].chunks.length;
        immutable num_chord_pt =  wing_state.wing_part_states[0].ctrl_chunks.length/num_span_chunks;
        double root_chord = wing.wing_parts[0].wing_root_chord; 
        size_t c_idx = 0;
        foreach(wp_idx, wing_part_state; wing_state.wing_part_states){
            foreach (span_idx; 0..num_span_chunks){
                wing_part_state.circulation_model.compute_dCl(wing_lift_surf,wing_part_state,wp_idx,span_idx);
                wing_part_state.chunks[span_idx].dC_L[] /= root_chord; //circulation is multiplied by root chord in compute_d_gamma_circulation, so Cl need to be devided by it
                dC_L[c_idx][] = wing_part_state.chunks[span_idx].dC_L[];
                c_idx++;
            }
        }

        wing_state.C_L = integrate_trapaziodal(dC_L, y);
    }

    InducedVelocities compute_wing_induced_vel_on_blade(immutable Chunk x, immutable Chunk y, immutable Chunk z){
        auto ind_vel = compute_wing_induced_vel(wing_lift_surf, x, y, z);

        return ind_vel;
    }
    
    Chunk inflow_at(immutable Chunk x, immutable Chunk y, immutable Chunk z, immutable Chunk x_e, double angle_of_attack){
        auto ind_vel = compute_wing_induced_vel(wing_lift_surf, x, y, z);

        immutable Chunk neg_z = ind_vel.v_z[];
        return neg_z;
    }

    Chunk inflow_at(immutable Vector!(4, Chunk) xyz) {
        auto ind_vel = compute_wing_induced_vel(wing_lift_surf, xyz[0], xyz[1], xyz[2]);

        immutable Chunk neg_z = ind_vel.v_z[];
        return neg_z;
    }
}