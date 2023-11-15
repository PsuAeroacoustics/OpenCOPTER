
import opencopter.aircraft;
import opencopter.atmosphere;
import opencopter.math;
import opencopter.memory;
import opencopter.vortexlattice;
import opencopter.config;
import opencopter.inflow;
import opencopter.wake;
import opencopter.airfoilmodels;
import opencopter.liftmodels;
import opencopter.bladeelement;


import std.algorithm;
import std.array;
import std.stdio;
import std.math;
import std.range: iota;
import std.conv;
import std.file : append;
import plt = matplotlibd.pyplot;

void main(){
    //wing properties
    double wing_span = 7.5;
    double le_sweep = 45.0; //degrees
    double te_sweep = 0.0; //degrees
    double wing_avg_chord = 1.0;
    double root_chord = 2.0;
    double tip_chord = 1.0;
    size_t chord_vortex_nodes = 5;
    size_t span_vortex_nodes = 20;
    //size_t num_rotors = 2;
    size_t num_wings = 1;
    size_t num_wing_parts = 2;
    double alpha = 7.0;
    Vec3 wing_origin = [0,0,0];
    //Vec3 rotor_origin_0 = [1,0,0];
    //Vec3 rotor_origin_1 = [-5,0,0];
    auto rotor_origin = [Vec3(1,0,0), Vec3(-5,0,0)];
    //auto rotor_origin = append(rotor_origin_0,rotor_origin_1);
    double freestream_velocity = 20.0;
    Chunk V_inf = 20.0;

    //rotor properties
    size_t iterations = 5400;
    size_t wake_history_length = 1*1024;
    size_t num_blade_elements = 45;

    size_t num_rotors = 2;
    size_t num_blades = 4;
    double R = 2.0;

    double theta_75 = 4.5*(PI/180);

    double density = 1.125;
    double omega = -109.12;
    double blade_AR = 16.5;
    double sos = 343;

    double rotor_aoa = -5.3*(PI/180);
    double theta_tw_1 = -8.0*(PI/180);

    double d_psi = 1.0;
    double dt = d_psi*(PI/180.0)/abs(omega);

    double shed_history_angle = 45.0;
    size_t shed_history = (shed_history_angle/d_psi).to!size_t;

    double azimuth_offset = 2.0*PI/num_blades;

    double[] r = generate_radius_points(num_blade_elements);
    size_t elements = r.length;
    double[] c = new double[elements];
    double[] alpha_0 = new double[elements];
    auto twist =  r.map!(x => (x-0.75)*theta_tw_1);
    double C_l_alpha = 2.0*PI;
    double[] sweep = new double[elements];

    auto atmo = Atmosphere(1.225,18.03e-6);

    auto airfoil = new ThinAirfoil(0.0401);

    AirfoilModel[] airfoils;

    airfoils[0] = airfoil;

    size_t[2][] extent = [[0,num_blade_elements-1]];
    //extent[0][0] = 0.0;
    //extent[1][0] = num_blade_elements-1;

    auto blade_airfoil = new BladeAirfoil(airfoils, extent);

    double rotor_avg_chord = R*sum(c)/c.length;

    //wing geometry
    double[] y = generate_spanwise_control_points(span_vortex_nodes);
    //writeln("y = ", y);
    
    double lamda = tip_chord/root_chord; //taper ratio;
    double AR = 2*wing_span/(tip_chord + root_chord);
    double delta_m = 4*(lamda- 1)/(AR*(lamda+ 1));
    auto chord_dist = y.map!(x => (1.0- (1.0-lamda)*2*x)).array;
    writeln("y_dist : ", y);    
    writeln("chord_dist : ", chord_dist);

    auto aircraft = Aircraft(num_rotors,num_wings);

    foreach(i;0..num_rotors){
        auto rotor = RotorGeometry(num_blades, rotor_origin[i], R, 0);
        foreach(j;0..num_blades){
            auto blade = BladeGeometry(num_blade_elements,azimuth_offset,rotor_avg_chord,blade_airfoil);
            rotor.blades[j] = blade;
        }
        rotor.solidity = num_blades*rotor.blades[0].average_chord/(PI*rotor.radius);
        aircraft.rotors[i] = rotor;
    }

    auto wing = WingGeometry(num_wing_parts, wing_origin, wing_span);

    aircraft.wings[0]=wing;
    
    size_t acutal_span_nodes = span_vortex_nodes%chunk_size == 0 ? span_vortex_nodes : span_vortex_nodes + (chunk_size - span_vortex_nodes%chunk_size);
    size_t span_chunks = acutal_span_nodes/chunk_size;

    foreach(i;0..num_wing_parts){
        auto wing_part = WingPartGeometry(span_vortex_nodes, chord_vortex_nodes, wing_origin, wing_avg_chord, root_chord, tip_chord, le_sweep, te_sweep, wing_span);
        //wing_part[i].set_geometry_array!"twist"(twist_dist);
        wing_part.set_geometry_array!"chord"(chord_dist);
        //wing_part[i].set_geometry_array!"C_l_alpha"(C_l_alpha_dist);
        //wing_part[i].set_geometry_array!"alpha_0"(alpha_0_dist);
        //wing_part[i].set_geometry_array!"sweep"(sweep_dist);

        aircraft.wings[0].wing_parts[i] = wing_part;
    }

    set_wing_ctrl_pt_geometry(aircraft.wings[0], span_vortex_nodes, chord_vortex_nodes);

    auto ac_state = AircraftState(num_rotors, 0, 0, num_wings, num_wing_parts, span_vortex_nodes, chord_vortex_nodes,aircraft);

    auto wing_lift_surface = WingLiftSurf(num_wing_parts);
    foreach(i;0..num_wing_parts){
        auto wing_part_lift_surf = WingPartLiftingSurf(span_vortex_nodes,chord_vortex_nodes);
        wing_lift_surface.wing_part_lift_surf[i] = wing_part_lift_surf;
    }
    

    set_wing_vortex_geometry(wing_lift_surface, aircraft.wings[0], span_chunks, chord_vortex_nodes);

    /*writeln("wing_vortex_geometry");
    foreach(wp_idx, wp_lift_surf; wing_lift_surface.wing_part_lift_surf){
        foreach(fl_idx, filament; wp_lift_surf.spanwise_filaments){
            foreach(ch_idx, chunk; filament.chunks){
                writeln("wp_idx = ", wp_idx, "fl_idx = ",fl_idx, "chunk_idx = ", ch_idx, "x = ", chunk.x);
                writeln("wp_idx = ", wp_idx, "fl_idx = ",fl_idx, "chunk_idx = ", ch_idx, "y = ", chunk.y);
                writeln("wp_idx = ", wp_idx, "fl_idx = ",fl_idx, "chunk_idx = ", ch_idx, "z = ", chunk.z);
            }        
        }
    }*/

    

    auto ac_input_state  = AircraftInputState(num_rotors,0,num_wings);
    ac_input_state.wing_inputs[0].freestream_velocity = freestream_velocity;

    foreach(i;0..num_rotors){
        ac_input_state.rotor_inputs[i].angle_of_attack = rotor_aoa;
        ac_input_state.rotor_inputs[i].angular_velocity = omega;
        ac_input_state.rotor_inputs[i].angular_accel = 0;
        ac_input_state.rotor_inputs[i].azimuth = 0;
        ac_input_state.rotor_inputs[i].freestream_velocity = freestream_velocity;
    }

    //inflows;

    long Mo = 4;
    long Me = 2;

    auto wake_history = WakeHistory(num_rotors,num_blades, wake_history_length,2 , num_blade_elements, shed_history);
    
    Inflow[] rotor_inflows = iota(0.0,num_rotors).map!((j){return cast(Inflow)(new HuangPetersInflow(Mo, Me, &aircraft.rotors[j.to!size_t], & ac_state.rotor_states[j.to!size_t], & ac_input_state.rotor_inputs[j.to!size_t], dt));}).array;

    Inflow[] wing_inflows = iota(0.0,num_wings).map!((j){return cast(Inflow)(new WingInflow(&aircraft.wings[j.to!size_t], &ac_state.wing_states[j.to!size_t], &ac_input_state.wing_inputs[j.to!size_t], &wing_lift_surface));}).array;
    
    auto all_inflows = rotor_inflows~wing_inflows;

    foreach(iteration; 0..iterations){
        if(iteration% 100 == 0){
            writeln("iterations : ", iteration);
        }
        foreach(r_idx, rotor_input; ac_input_state.rotor_inputs){
            rotor_input.azimuth += rotor_input.angular_velocity*dt + rotor_input.angular_accel*dt*dt;

            if(rotor_input.azimuth > 2.0*PI){
                rotor_input.azimuth = rotor_input.azimuth % (2*PI);
            } 
        }

        step(ac_state, aircraft, ac_input_state, all_inflows, wake_history, atmo, iteration, dt);

        if(iteration%100 == 0){
            writeln("C_T_0 = ", ac_state.rotor_states[0].C_T, "\tC_T_1= ", ac_state.rotor_states[1].C_T);
        }
    }


    //Wing alone simulation (tested already)
    foreach(wp_idx, wing_part_state; ac_state.wing_states[0].wing_part_states){
        foreach(sp_chunk; 0..chord_vortex_nodes*span_chunks){
            foreach(c1;0..chunk_size){
                wing_part_state.ctrl_chunks[sp_chunk].ctrl_pt_aoa[c1] = alpha*PI/180;
            }
        }
    }

        
    foreach(wp_idx, wing_part_state; ac_state.wing_states[0].wing_part_states){
        foreach(chord_idx; 0..chord_vortex_nodes){
            foreach(ch_idx;0..span_chunks){
                wing_part_state.circulation_model.compute_d_gamma_coefficients(wing_lift_surface, wing_part_state, wp_idx, ch_idx, chord_idx, V_inf[]);
                writeln("wing_part = ", wp_idx, "\tchord_idx = ", chord_idx, "\tspan_chunk = ch_idx ", ch_idx,"\tA_lk = ",wing_lift_surface.wing_part_lift_surf[wp_idx].spanwise_filaments[chord_idx].chunks[ch_idx].A_kl);
                //wing_part_state.circulation_model.compute_bound_circulation(wing_lift_surface, aircraft.wings[0].wing_parts[wp_idx], wp_idx, ch_idx, chord_idx);
                //writeln("wing_part = ", wp_idx, "\tchord_idx = ", chord_idx, "\tspan_chunk = ch_idx ", ch_idx,"\tA_lk = ",wing_lift_surface.wing_part_lift_surf[wp_idx].spanwise_filaments[chord_idx].chunks[ch_idx].gamma);
            //double[] A_kl_print = wing_lift_surface.spanwise_filaments[chord_idx].chunks[ch_idx].A_kl[].array;
            //writeln(A_kl_print);
                }
            //wing_part_state.circulation_model.compute_dCl(wing_lift_surface, wing_part_state, wp_idx, ch_idx, V_inf[]);
        }
    }
    writeln("\n");
    foreach(wp_idx, wing_part_state; ac_state.wing_states[0].wing_part_states){
        foreach(chord_idx; 0..chord_vortex_nodes){
            foreach(ch_idx;0..span_chunks){
                //wing_part_state.circulation_model.compute_d_gamma_coefficients(wing_lift_surface, wing_part_state, wp_idx, ch_idx, chord_idx, V_inf[]);
                //writeln("wing_part = ", wp_idx, "\tchord_idx = ", chord_idx, "\tspan_chunk = ch_idx ", ch_idx,"\tA_lk = ",wing_lift_surface.wing_part_lift_surf[wp_idx].spanwise_filaments[chord_idx].chunks[ch_idx].A_kl);
                wing_part_state.circulation_model.compute_bound_circulation(wing_lift_surface, aircraft.wings[0].wing_parts[wp_idx], wp_idx, ch_idx, chord_idx);
                //writeln("wing_part = ", wp_idx, "\tchord_idx = ", chord_idx, "\tspan_chunk = ch_idx ", ch_idx,"\tgamma = ",wing_lift_surface.wing_part_lift_surf[wp_idx].spanwise_filaments[chord_idx].chunks[ch_idx].gamma);
            //double[] A_kl_print = wing_lift_surface.spanwise_filaments[chord_idx].chunks[ch_idx].A_kl[].array;
            //writeln(A_kl_print);
                }
            //wing_part_state.circulation_model.compute_dCl(wing_lift_surface, wing_part_state, wp_idx, ch_idx, V_inf[]);
        }
    }

    foreach(wp_idx, wing_part_state; ac_state.wing_states[0].wing_part_states){
        foreach(ch_idx;0..span_chunks){
            wing_part_state.circulation_model.compute_dCl(wing_lift_surface, wing_part_state, wp_idx, ch_idx, V_inf[]);
        }
    }

    //writeln("chord_idx = 1\t right_wing,\tgamma = ",wing_lift_surface.wing_part_lift_surf[1].spanwise_filaments[0].chunks[1].gamma);
    writeln("chord_idx = 1\t right_wing,\td_Cl = ",ac_state.wing_states[0].wing_part_states[0].chunks[0].dC_L);


    double[] sectional_Cl = new double[span_chunks*chunk_size];
    writeln("length of sectional Cl array= ",sectional_Cl[].length);
    foreach(span_idx; 0..span_chunks){
        foreach(c1; 0..chunk_size){
            sectional_Cl[span_idx*chunk_size..(span_idx+1)*chunk_size] = ac_state.wing_states[0].wing_part_states[1].chunks[span_idx].dC_L;
        } 
    }


    Chunk x_pt = 5.0;
    Chunk z_pt = 2.0;
    double[8] y_pt = [-4.0,-3.0,-2.0,-1.0,0.0,1.0,2.0,3.0];
    double[] v_x = new double[y_pt.length];
    double[] v_y = new double[y_pt.length];
    double[] v_z = new double[y_pt.length];

    foreach(ch1; 0..y_pt.length/chunk_size){
        auto ind_vel = compute_wing_induced_vel(wing_lift_surface, x_pt, y_pt, z_pt);
        v_x[ch1*chunk_size..ch1*chunk_size + chunk_size] = ind_vel.v_x;
        v_y[ch1*chunk_size..ch1*chunk_size + chunk_size] = ind_vel.v_y;
        v_z[ch1*chunk_size..ch1*chunk_size + chunk_size] = ind_vel.v_z;
    }

    writeln("v_x = ", v_x);
    writeln("v_y = ", v_y);
    writeln("v_z = ", v_z);


    plt.figure();
    plt.plot(y.map!(x => x*2),sectional_Cl[]);
    plt.xlabel("2y/b");
    plt.ylabel("sectional Cl");
    plt.grid();
    plt.savefig("Cl_left_wing.png");

    /*plt.figure();
    plt.plot(y[],sectional_Cl[][1]);
    plt.xlabel("2y/b");
    plt.ylabel("sectional Cl");
    plt.grid();
    plt.savefig("Cl_right_wing.png");*/

    plt.figure();
    plt.plot(y_pt, v_x[], ["label": "v_x"]);
    plt.plot(y_pt, v_y[], ["label": "v_y"]);
    plt.plot(y_pt, v_z[], ["label": "v_z"]);
    plt.xlabel("y");
    plt.ylabel("Induced velocity");
    plt.grid();
    plt.savefig("Wing induced velocity");

}







