
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
    double wing_span = 7.5;  //half wing span
    double le_sweep = 45.0; //degrees
    double te_sweep = 0.0; //degrees
    double wing_avg_chord = 1.5;
    double root_chord = 2;
    double tip_chord = 1;
    size_t chord_vortex_nodes = 4;
    size_t span_vortex_nodes = 8;
    //size_t num_rotors = 2;
    size_t num_wings = 1;
    size_t num_wing_parts = 2;
    double wing_alpha = 7.0;
    Vec3 wing_origin = Vec3(0.0,0.0,0.0);//Vec3((0.102)*0.8255,(-0.215)*0.8255 - 0.3048,-0.285*0.8255);
    //Vec3 rotor_origin_0 = [1,0,0];
    //Vec3 rotor_origin_1 = [-5,0,0];
    auto rotor_origin = [Vec3(0.0,0.0,0.0)];
    //auto rotor_origin = append(rotor_origin_0,rotor_origin_1);
    double freestream_velocity = 25.7;
    Chunk V_inf = 25.7;

    //rotor properties
    size_t iterations = 3600;
    size_t wake_history_length = 1*1024;
    size_t num_blade_elements = 48;

    size_t num_rotors = 0;
    size_t num_blades = 0;
    double R = 0.8255;

    double theta_75 = 9.0*(PI/180);

    double density = 1.125;
    double omega = 207.345;
    double blade_AR = 12.99;
    double sos = 343;

    double rotor_aoa = 6.0*(PI/180);
    double theta_tw_1 = -13.0*(PI/180);

    double d_psi = 1.0;
    double dt = d_psi*(PI/180.0)/abs(omega);

    double shed_history_angle = 45.0;
    size_t shed_history = (shed_history_angle/d_psi).to!size_t;

    double azimuth_offset = 2.0*PI/num_blades;

    double[] r = generate_radius_points(num_blade_elements);
    size_t elements = r.length;
    double[] c = r.map!(x => 1/blade_AR).array;
    double[] alpha_0 = new double[elements];
    double[] twist =  r.map!(x => (x-0.75)*theta_tw_1).array;
    double[] C_l_alpha = r.map!(x => 2.0*PI.to!double).array;
    double[] sweep = r.map!(x => 0.0).array;

    auto atmo = Atmosphere(1.225,18.03e-6);

    auto airfoil = new ThinAirfoil(0.0401);

    AirfoilModel[1] airfoils = [airfoil];

    //airfoils[0] = airfoil;

    size_t[2][] extent = [[0,num_blade_elements-1]];
    //extent[0][0] = 0.0;
    //extent[1][0] = num_blade_elements-1;

    auto blade_airfoil = new BladeAirfoil(airfoils, extent);

    double rotor_avg_chord = R*sum(c)/c.length;

    size_t iter_per_revs = (360/d_psi).to!size_t;
    size_t num_revs = iterations/iter_per_revs;

    //wing geometry
    double[] y = generate_spanwise_control_points(span_vortex_nodes);
    //writeln("y = ", y);
    
    double lamda = tip_chord/root_chord; //taper ratio;
    double AR = 2*wing_span/(tip_chord + root_chord);
    double delta_m = 4*(lamda- 1)/(AR*(lamda+ 1));
    auto chord_dist = y.map!(x => (1.0- (1.0-lamda)*2*x)).array;
    double camber = 0.0;

    /*foreach(chord_idx; 0..chord_vortex_nodes){
        camber[chord_idx][] = y.map!(x => 0.0).array;
    }*/
    //writeln("y_dist : ", y);    
    //writeln("chord_dist : ", chord_dist);

    // Output data
    double[][] rotor_CT = allocate_dense(iter_per_revs,num_blade_elements);
    double[][] wing_part_dCl = allocate_dense(num_wing_parts,span_vortex_nodes);
    double[][] wing_dCl = allocate_dense(iter_per_revs,2*span_vortex_nodes);
    double[] wing_dCl_root = new double[iter_per_revs];
    double[] wing_dCl_mid = new double[iter_per_revs];
    double[] wing_dCl_tip = new double[iter_per_revs];
    double[] rotor_CT_mid = new double[iter_per_revs];
    double[][] ind_vel_at_wing = new double[][2*span_vortex_nodes];
    double[][] ind_vel_at_rotor = new double[][4];

    auto aircraft = Aircraft(num_rotors,num_wings);

    /*foreach(i;0..num_rotors){
        auto rotor = RotorGeometry(num_blades, rotor_origin[i], R, 0);
        foreach(j;0..num_blades){
            auto blade = BladeGeometry(num_blade_elements,j*azimuth_offset,rotor_avg_chord,blade_airfoil);
            blade.set_geometry_array!"r"(r);
            blade.set_geometry_array!"twist"(twist);
            blade.set_geometry_array!"alpha_0"(alpha_0);
            blade.set_geometry_array!"C_l_alpha"(C_l_alpha);
            blade.set_geometry_array!"chord"(c);
            blade.set_geometry_array!"sweep"(sweep);

            rotor.blades[j] = blade;
        }
        rotor.solidity = num_blades*rotor.blades[0].average_chord/(PI*rotor.radius);
        aircraft.rotors[i] = rotor;
    }*/

    //writeln("blade in each rotor = ", aircraft.rotors[0].blades.length);
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

    set_wing_ctrl_pt_geometry(aircraft.wings[0], span_vortex_nodes, chord_vortex_nodes, camber);

    

    auto wing_lift_surface = WingLiftSurf(num_wing_parts);
    foreach(i;0..num_wing_parts){
        auto wing_part_lift_surf = WingPartLiftingSurf(span_vortex_nodes,chord_vortex_nodes);
        wing_lift_surface.wing_part_lift_surf[i] = wing_part_lift_surf;
    }
    

    set_wing_vortex_geometry(wing_lift_surface, aircraft.wings[0], span_chunks, chord_vortex_nodes);

    set_circulation_to_zero(wing_lift_surface);
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

    writeln("instantiating ac_state");

    auto ac_state = AircraftState(num_rotors, num_blades, num_blade_elements, num_wings, num_wing_parts, span_vortex_nodes, chord_vortex_nodes,aircraft);
    writeln("instantiating ac_state : done");

    auto ac_input_state  = AircraftInputState(num_rotors,num_blades,num_wings);
    ac_input_state.wing_inputs[0].freestream_velocity = freestream_velocity;

    /*foreach(r_idx;0..num_rotors){
        ac_input_state.rotor_inputs[r_idx].angle_of_attack = rotor_aoa;
        ac_input_state.rotor_inputs[r_idx].cos_aoa = cos(ac_input_state.rotor_inputs[r_idx].angle_of_attack);
		ac_input_state.rotor_inputs[r_idx].sin_aoa = sin(ac_input_state.rotor_inputs[r_idx].angle_of_attack);
        ac_input_state.rotor_inputs[r_idx].angular_velocity = omega;
        ac_input_state.rotor_inputs[r_idx].angular_accel = 0;
        ac_input_state.rotor_inputs[r_idx].azimuth = 0;
        ac_input_state.rotor_inputs[r_idx].freestream_velocity = freestream_velocity;
        ac_state.rotor_states[r_idx].advance_ratio = ac_input_state.rotor_inputs[r_idx].freestream_velocity * ac_input_state.rotor_inputs[r_idx].cos_aoa/abs(ac_input_state.rotor_inputs[r_idx].angular_velocity*aircraft.rotors[r_idx].radius);
        ac_state.rotor_states[r_idx].axial_advance_ratio = ac_input_state.rotor_inputs[r_idx].freestream_velocity * ac_input_state.rotor_inputs[r_idx].sin_aoa/abs(ac_input_state.rotor_inputs[r_idx].angular_velocity*aircraft.rotors[r_idx].radius);
        //writeln("at start of the simulation, ", "adv_ratio=", ac_state.rotor_states[r_idx].advance_ratio, "axial_adv_ratio=", ac_state.rotor_states[r_idx].axial_advance_ratio);
        foreach(b_idx;0..num_blades){
            ac_input_state.rotor_inputs[r_idx].r_0[b_idx] = 0.1*aircraft.rotors[r_idx].blades[b_idx].average_chord;
            ac_input_state.rotor_inputs[r_idx].blade_flapping[b_idx] = 0.0;
            ac_input_state.rotor_inputs[r_idx].blade_flapping_rate[b_idx] = 0.0;
            ac_input_state.rotor_inputs[r_idx].blade_pitches[b_idx] = theta_75;
        }
    }*/

    ac_input_state.wing_inputs[0].angle_of_attack = wing_alpha;
    //ac_input_state.wing_inputs[0].freestream_velocity = freestream_velocity;


    
    //inflows;

    /*long Mo = 4;
    long Me = 2;

    auto wake_history = WakeHistory(num_rotors,num_blades, wake_history_length, 2 , num_blade_elements, shed_history);
    
    Inflow[] rotor_inflows = iota(0.0,num_rotors).map!((j){return cast(Inflow)(new HuangPetersInflow(Mo, Me, &aircraft.rotors[j.to!size_t], & ac_state.rotor_states[j.to!size_t], & ac_input_state.rotor_inputs[j.to!size_t], dt));}).array;

    Inflow[] wing_inflows = iota(0.0,num_wings).map!((j){return cast(Inflow)(new WingInflow(&aircraft.wings[j.to!size_t], &ac_state.wing_states[j.to!size_t], &ac_input_state.wing_inputs[j.to!size_t], &wing_lift_surface));}).array;
    
    auto all_inflows = rotor_inflows~wing_inflows;
    size_t wrt_iter = 0;

    //File f = File("data.dat", "w+");
    foreach(iteration; 0..iterations){
        //if(iteration% 100 == 0){
           
        //}
        foreach(r_idx, rotor_input; ac_input_state.rotor_inputs){
            rotor_input.azimuth += rotor_input.angular_velocity*dt + rotor_input.angular_accel*dt*dt;

            if(rotor_input.azimuth > 2.0*PI){
                rotor_input.azimuth = rotor_input.azimuth % (2*PI);
            } 
        }
        //writeln("going to step function, iteration: ", iteration);
        step(ac_state, aircraft, ac_input_state, all_inflows, wake_history, atmo, iteration, dt);
        if(iteration%iter_per_revs == 0){
            writeln("iterations : ", iteration, "\tC_T_0 = ", ac_state.rotor_states[0].C_T );
        }

        if((iteration- iteration%iter_per_revs)/iter_per_revs == num_revs-1){
            foreach(ch1; 0..ac_state.rotor_states[0].blade_states[0].chunks.length){
                rotor_CT[wrt_iter][ch1*chunk_size..ch1*chunk_size + chunk_size] = ac_state.rotor_states[0].blade_states[0].chunks[ch1].dC_T;
            }       

        foreach(wp_idx, wing_part; ac_state.wing_states[0].wing_part_states){
            foreach(ch1; 0..wing_part.chunks.length){
                wing_part_dCl[wp_idx][ch1*chunk_size..ch1*chunk_size+ chunk_size] = wing_part.chunks[ch1].dC_L;
            }
        }
            wing_dCl_mid[wrt_iter] = wing_part_dCl[1][0];
            wing_dCl_root[wrt_iter] = wing_part_dCl[0][7];
            wing_dCl_tip[wrt_iter] = wing_part_dCl[1][7];
            rotor_CT_mid[wrt_iter] = rotor_CT[wrt_iter][24];
            //writeln("iteration = ", wrt_iter, "wing_dCl = ", wing_dCl[wrt_iter][]);
            wrt_iter++;
        }      
        
    }

    double[] azimuth = new double[360];
    foreach(i; 0..azimuth.length){
        azimuth[i] = i;
    }

    plt.figure();
    plt.plot(r,rotor_CT[359][0..$]);
    plt.xlabel("r");
    plt.ylabel("C_l");
    plt.grid();
    plt.savefig("CT_blade_last_iteration.png");

    plt.figure();
    plt.plot(azimuth,rotor_CT_mid);
    plt.xlabel("azimuth");
    plt.ylabel("C_l");
    plt.grid();
    plt.savefig("CT_blade_mid.png");

    plt.figure();
    plt.plot(azimuth,wing_dCl_root[]);
    plt.xlabel("azimuth");
    plt.ylabel("C_l");
    plt.grid();
    plt.savefig("Cl_wing_root.png");

    plt.figure();
    plt.plot(azimuth,wing_dCl_mid[]);
    plt.xlabel("azimuth");
    plt.ylabel("C_l");
    plt.grid();
    plt.savefig("Cl_wing_mid.png");

    plt.figure();
    plt.plot(azimuth,wing_dCl_tip[]);
    plt.xlabel("azimuth");
    plt.ylabel("C_l");
    plt.grid();
    plt.savefig("Cl_wing_tip.png");
    */

    

    


    //Wing alone simulation (tested already)
    foreach(wp_idx, wing_part_state; ac_state.wing_states[0].wing_part_states){
        foreach(sp_chunk; 0..chord_vortex_nodes*span_chunks){
            foreach(c1;0..chunk_size){
                wing_part_state.ctrl_chunks[sp_chunk].ctrl_pt_aoa[c1] = wing_alpha*PI/180;
                wing_part_state.ctrl_chunks[sp_chunk].ctrl_pt_ut[c1] = freestream_velocity;
            }
        }
    }

        
    foreach(wp_idx, wing_part_state; ac_state.wing_states[0].wing_part_states){
        foreach(chord_idx; 0..chord_vortex_nodes){
            foreach(ch_idx;0..span_chunks){
                wing_part_state.circulation_model.compute_d_gamma_coefficients(wing_lift_surface, wing_part_state, wp_idx, ch_idx, chord_idx, V_inf[]);
                //writeln("wing_part = ", wp_idx, "\tchord_idx = ", chord_idx, "\tspan_chunk = ch_idx ", ch_idx,"\tA_lk = ",wing_lift_surface.wing_part_lift_surf[wp_idx].spanwise_filaments[chord_idx].chunks[ch_idx].A_kl);
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
                writeln("wing_part = ", wp_idx, "\tchord_idx = ", chord_idx, "\tspan_chunk = ch_idx ", ch_idx,"\tgamma = ",wing_lift_surface.wing_part_lift_surf[wp_idx].spanwise_filaments[chord_idx].chunks[ch_idx].gamma);
            //double[] A_kl_print = wing_lift_surface.spanwise_filaments[chord_idx].chunks[ch_idx].A_kl[].array;
            //writeln(A_kl_print);
                }
            //wing_part_state.circulation_model.compute_dCl(wing_lift_surface, wing_part_state, wp_idx, ch_idx, V_inf[]);
        }
    }

    foreach(wp_idx, wing_part_state; ac_state.wing_states[0].wing_part_states){
        foreach(ch_idx;0..span_chunks){
            wing_part_state.circulation_model.compute_dCl(wing_lift_surface, wing_part_state, wp_idx, ch_idx);
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


    /*plt.figure();
    plt.plot(y.map!(x => x*2),sectional_Cl[]);
    plt.xlabel("2y/b");
    plt.ylabel("sectional Cl");
    plt.grid();
    plt.savefig("Cl_left_wing.png");*/

    /*plt.figure();
    plt.plot(y[],sectional_Cl[][1]);
    plt.xlabel("2y/b");
    plt.ylabel("sectional Cl");
    plt.grid();
    plt.savefig("Cl_right_wing.png");

    plt.figure();
    plt.plot(y_pt, v_x[], ["label": "v_x"]);
    plt.plot(y_pt, v_y[], ["label": "v_y"]);
    plt.plot(y_pt, v_z[], ["label": "v_z"]);
    plt.xlabel("y");
    plt.ylabel("Induced velocity");
    plt.grid();
    plt.savefig("Wing induced velocity");*/

}







