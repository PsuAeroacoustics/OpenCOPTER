
import sys
import os

# Tell pythn where it will find the libopencopter module
sys.path.append(f'{os.path.dirname(os.path.realpath(__file__))}/../')

from libopencopter import *
import numpy as np
import math
import time
import scipy.io
import matplotlib.pyplot as plt

if __name__ == "__main__":
    
    # Wing geometric properties
    wing_span = 0.3048 # this is half wing span (b/2)
    le_sweep = 0.0
    te_sweep = 0.0
    wing_avg_chord = 0.3048
    root_chord = 0.3048
    tip_chord = 0.3048
    
    chord_vortex_nodes = 4
    span_vortex_nodes = 8
    num_wings = 1
    num_wing_parts = [2]
    wing_alpha = 0.0
    wing_origin = [Vec3([(0.102)*0.8255, ((-0.215)*0.8255 - 0.3048), -0.285*0.8255])]
    
    rotor_origin = [Vec3([0.0, 0.0, 0.0])]
    freestream_velocity = 25.7
    
    # Rotor geometric properties
    iterations = 3600
    wake_history_length = [1*1024]
    num_blade_elements = 48
    
    num_rotors = 1
    num_blades = [4]
    rotor_R = 0.8255
    
    theta_75 = 9.0*(math.pi/180.0)
    
    density = 1.125
    omega = 207.345
    blade_AR = 12.99
    sos = 343
    
    rotor_aoa = 6.0*(math.pi/180)
    theta_tw_1 = -13.0*(math.pi/180)
    
    d_psi = 1.0 # degrees
    dt = d_psi*(math.pi/180.0)/abs(omega)
    
    shed_history_angle = 45.0
    shed_history = round(shed_history_angle/d_psi)
    
    d_azimuth = 2.0*math.pi/num_blades[0]
    
    # Create the spanwise distributions for our blades
	# We will later apply these to opencopters internal
	# data layout. r and c arrays are normalized by the
	# rotor radius
    r = generate_radius_points(num_blade_elements)
    elements = len(r)
    print("requested_elements: ", num_blade_elements, " actual elements: ", elements)
    c = (1.0/blade_AR)*np.ones(elements)
    alpha_0 = np.zeros(elements)
    twist = [(_r - 0.75)*theta_tw_1 for _r in r]
    C_l_alpha = 2.0*math.pi*np.ones(elements)
    sweep = np.zeros(elements)

    atmo = Atmosphere(1.125, 18.03e-6)
    
    extent = [[0,num_blade_elements-1]]
    airfoil = [ThinAirfoil(0.0401)]
    
    blade_airfoil = BladeAirfoil(airfoil,extent)
    
    avg_rotor_chord = rotor_R*sum(c)/len(c)
    
    iter_per_revs = (360/d_psi)
    num_revs = iterations/iter_per_revs
    
    wing_y = generate_spanwise_control_points(span_vortex_nodes)
    
    lamda = tip_chord/root_chord
    AR = 2*wing_span/(tip_chord + root_chord)
    delta_m = 4*(lamda - 1)/(AR*(lamda + 1))
    chord_distribution = [(1.0 - (1.0-lamda)*2.0*y)*root_chord for y in wing_y]
    camber = 0.0
    wing_twist = np.zeros(span_vortex_nodes)
    wing_sweep = np.zeros(span_vortex_nodes) # change it to a formula based on wing parameters
    
    aircraft = Aircraft(num_rotors, num_wings)
    
    def build_blade(b_idx):
        
        blade = BladeGeometry(
            num_elements = elements,
            azimuth_offset = b_idx*d_azimuth,
            average_chord = rotor_R*np.sum(c)/len(c),
            airfoil = blade_airfoil
        )
        
        set_r(blade, r)
        set_twist(blade, twist)
        set_alpha_0(blade, alpha_0)
        set_C_l_alpha(blade, C_l_alpha)
        set_chord(blade, c)
        set_sweep(blade, sweep)
        
        return blade
    
    def build_rotor(r_idx):
        rotor = RotorGeometry(
            num_blades = num_blades[r_idx],
            origin = rotor_origin[r_idx],
            radius = rotor_R,
            solidity = 0
        )
        
        rotor.blades = [build_blade(b_idx) for b_idx in range(num_blades[r_idx])]
        rotor.solidity = num_blades[r_idx]*rotor.blades[0].average_chord/(math.pi*rotor.radius)
        return rotor

    aircraft.rotors = [build_rotor(r_idx) for r_idx in range(num_rotors)]
    
    
    def build_wing_part(wp_idx):
        
        wing_part = WingPartGeometry(
            span_elements = span_vortex_nodes,
            chordwise_nodes = chord_vortex_nodes,
            wing_root_origin = wing_origin[0],
            average_chord = wing_avg_chord,
            wing_root_chord = root_chord,
            wing_tip_chord = tip_chord,
            le_sweep_angle = le_sweep,
            te_sweep_angle = te_sweep,
            wing_span = wing_span
        )
        
        set_wing_chord(wing_part,chord_distribution)
        set_wing_twist(wing_part,wing_twist)
        set_wing_sweep(wing_part,wing_sweep)
        
        return wing_part
    
    def build_wing(w_idx):
        wing = WingGeometry(
            num_parts = num_wing_parts[w_idx],
            origin = wing_origin[w_idx],
            wing_span = wing_span
        )

        wing.wing_parts = [build_wing_part(wp_idx) for wp_idx in range(num_wing_parts[w_idx])]
        
        return wing
    
    aircraft.wings= [build_wing(w_idx) for w_idx in range(num_wings)]
    
    set_wing_ctrl_pt_geometry(aircraft.wings[0], span_vortex_nodes, chord_vortex_nodes, camber)
    
    ac_state = AircraftState(num_rotors, num_blades, elements, num_wings, num_wing_parts, span_vortex_nodes, chord_vortex_nodes, aircraft)
    
    print("Freestream vel: ", freestream_velocity, "m/s")
    
    wake_history = WakeHistory(num_rotors, num_blades, wake_history_length, 2, elements, shed_history)
    
    wing_lift_surface = WingLiftSurf(num_wing_parts[0])
    for wp_idx in range(num_wing_parts[0]):
        wing_part_lift_surf = WingPartLiftingSurf(span_vortex_nodes, chord_vortex_nodes)
        wing_lift_surface.wing_part_lift_surf[wp_idx] = wing_part_lift_surf
    
    set_wing_vortex_geometry(wing_lift_surface, aircraft.wings[0], span_vortex_nodes, chord_vortex_nodes)
    
    ac_input_state = AircraftInputState(num_rotors, num_blades, num_wings)

    ac_input_state.rotor_inputs[0].angle_of_attack = rotor_aoa
    ac_input_state.rotor_inputs[0].angular_velocity = omega
    ac_input_state.rotor_inputs[0].angular_accel = 0
    ac_input_state.rotor_inputs[0].azimuth = 0
    ac_input_state.rotor_inputs[0].freestream_velocity = freestream_velocity
    for b_idx in range(num_blades[0]):
        ac_input_state.rotor_inputs[0].r_0[b_idx] = 0.1*aircraft.rotors[0].blades[b_idx].average_chord
        ac_input_state.rotor_inputs[0].blade_flapping[b_idx] = 0
        ac_input_state.rotor_inputs[0].blade_flapping_rate[b_idx] = 0
        ac_input_state.rotor_inputs[0].blade_pitches[b_idx] = theta_75
        
    ac_input_state.wing_inputs[0].angle_of_attack = wing_alpha
    ac_input_state.wing_inputs[0].freestream_velocity = freestream_velocity
    
    r_0 = ac_input_state.rotor_inputs[0].r_0[b_idx]
    print("r_0: ", ac_input_state.rotor_inputs[0].r_0)
    
    rotor_inflows = [HuangPeters(4, 2, aircraft.rotors[r_idx], ac_state.rotor_states[r_idx], ac_input_state.rotor_inputs[r_idx], dt) for r_idx in range(num_rotors)]
    wing_inflows = [WingInflow(aircraft.wings[w_idx], ac_state.wing_states[w_idx], ac_input_state.wing_inputs[w_idx], wing_lift_surface) for w_idx in range(num_wings)]
    
    inflows = np.append(rotor_inflows, wing_inflows)
    
    start_time = time.monotonic()
    for iteration in range(iterations):
        if iteration % 100 == 0:
            now = time.monotonic()
            print(now - start_time, ": iteration: ", iteration)
            start_time = now
        
        for r_idx in range(ac_input_state.rotor_inputs.length()):
            ac_input_state.rotor_inputs[r_idx].azimuth += ac_input_state.rotor_inputs[r_idx].angular_velocity*dt + ac_input_state.rotor_inputs[r_idx].angular_accel*dt*dt

            if ac_input_state.rotor_inputs[r_idx].azimuth > 2.0*math.pi:
                ac_input_state.rotor_inputs[r_idx].azimuth = math.fmod(ac_input_state.rotor_inputs[r_idx].azimuth, 2.0*math.pi)

        step(ac_state, aircraft, ac_input_state, inflows, wake_history, atmo, iteration, dt)
        
        if( iteration%iter_per_revs == 0):
            print("iterations : ", iteration, "\tC_T_0 = ", ac_staate.rotor_states[0].C_T)
        
        if iteration > (iterations - 360):
            write_rotor_vtu("rotor", iteration, 0, vtk_rotor, ac_state.rotor_states[0], ac_input_state.rotor_inputs[0])
            write_wake_vtu("wake", iteration, vtk_wake, wake_history.history[0])
            
        result_dictionary = {}
        spanwise_chunks = span_vortex_nodes/chunk_size()
        if(iteration == iterations-1):
            for wp_idx in range(num_wing_parts):
                wing_part_C_l = np.zeros(span_vortex_nodes)
                wing_part_gamma = np.zeros(span_vortex_nodes, chord_vortex_nodes)
                wing_part_A_kl = np.zeros(span_vortex_nodes, chord_vortex_nodes)
            
                for i in range(spanwise_chunks):
                    wing_part_C_l[i: (i+1)*chunk_size()] = ac_state.wings[0].wing_part_states[wp_idx].dC_L
                    for j in range(chord_vortex_nodes):
                        wing_part_gamma[i: (i+1)*chunk_size(), j] = wing_lift_surface.wing_part_lift_surf[wp_idx].spanwise_filaments[j].chunks[i].gamma
                        wing_part_A_kl[i: (i+1)*chunk_size(), j] = wing_lift_surface.wing_part_lift_surf[wp_idx].spanwise_filaments[j].chunks[i].A_kl
            
                result_dictionary[f'wing_part_{wp_idx}_C_l'] = wing_part_C_l
                result_dictionary[f'wing_part_{wp_idx}_gamma'] = wing_part_gamma
                result_dictionary[f'wing_part_{wp_idx}_A_kl'] - wing_part_A_kl
        
            for r_idx in range(num_rotors):
                actual_wake_history = wake_history_length[r_idx] if wake_history_length[r_idx]%chunk_size() == 0 else wake_history_length[r_idx] + (chunk_size() - wake_history_length[r_idx]%chunk_size())
                wake_trajectories = np.zeros((num_blades, 3, actual_wake_history))

                for b_idx in range(num_blades):
                    wake_trajectories[b_idx, 0, :] = get_wake_x_component(rotor_wake_history.history[0].rotor_wakes[r_idx].tip_vortices[b_idx])
                    wake_trajectories[b_idx, 1, :] = get_wake_y_component(rotor_wake_history.history[0].rotor_wakes[r_idx].tip_vortices[b_idx])
                    wake_trajectories[b_idx, 2, :] = get_wake_z_component(rotor_wake_history.history[0].rotor_wakes[r_idx].tip_vortices[b_idx])

                result_dictionary[f'wake_{r_idx}_trajectory'] = wake_trajectories

            result_dictionary["rotor_c_t"] = rotorcraft_thrusts
            result_dictionary["rotor_collectives"] = [rotorcraft_input_state.rotor_inputs[r_idx].blade_pitches[0] for r_idx in range(num_rotors)]
            result_dictionary["rotor_chis"] = [rotorcraft_inflows[r_idx].wake_skew() for r_idx in range(num_rotors)]
        
        if((iterations - iteration)<= iter_per_revs):
            blade_C_T = np.zeros(num_blades,iter_per_revs, elements)
            wing_C_L = np.zeros(num_wing_parts,iter_per_revs, span_vortex_nodes)
            
            for b_idx in range(num_blades):
                for ch_idx, chunk in enumerate(ac_state.rotor_states[0].blade_states[b_idx].chunks):
                    blade_C_T[b_idx, iteration, ch_idx: (ch_idx+1)*chunk_size()] = chunk.dC_T
            
            for wp_idx in range(num_wing_parts):
                for ch_idx, chunk in enumerate(ac_state.wing_states[0].wing_part_states[wp_idx].chunks):
                    wing_C_L[ wp_idx, iteration, ch_idx: (ch_idx+1)*chunk_size()] = chunk.dC_L
        
    scipy.io.savemat(f"/results.mat", result_dictionary)
                    
    print("rotor 0 C_T: ", ac_state.rotor_states[0].C_T)

    max_dim_1 = -math.inf
    min_dim_1 = math.inf

    max_dim_2 = -math.inf
    min_dim_2 = math.inf

    # Iterating directly on the collections like this is read only
    for r_idx, rotor in enumerate(ac_state.rotor_states):
        for b_idx, blade in enumerate(rotor.blade_states):
            x = get_wake_x_component(wake_history.history[0].rotor_wakes[r_idx].tip_vortices[b_idx])
            z = get_wake_z_component(wake_history.history[0].rotor_wakes[r_idx].tip_vortices[b_idx])

            max_dim_2 = max(max_dim_2, max(z))
            min_dim_2 = min(min_dim_2, min(z))

            max_dim_1 = max(max_dim_1, max(x))
            min_dim_1 = min(min_dim_1, min(x))
            plt.plot(x, z, linewidth=0.5)

            x_r = math.cos(blade.azimuth)
            x_b = aircraft.rotors[r_idx].origin[0] + x_r*math.cos(ac_input_state.rotor_inputs[r_idx].angle_of_attack)
            z_b = aircraft.rotors[r_idx].origin[2] - x_r*math.sin(ac_input_state.rotor_inputs[r_idx].angle_of_attack)

            plt.plot([aircraft.rotors[r_idx].origin[0], x_b], [aircraft.rotors[r_idx].origin[2], z_b], "k-", linewidth=1.5)
    
    plt.axis("square")

    span_1 = max_dim_1 - min_dim_1
    span_2 = max_dim_2 - min_dim_2
    plt.xlim(left=min_dim_1 - .1*span_1, right= max_dim_1 + .1*span_1)
    plt.ylim(bottom=min_dim_2 - .9*span_2, top=max_dim_2 + .9*span_2)
    plt.show()
    