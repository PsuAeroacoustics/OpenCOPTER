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

    wing_span = 7.5
    le_sweep = 45.0*math.pi/180
    te_sweep = 0.0
    span_vortex_nodes = 24
    chord_vortex_nodes = 5
    alpha = 7*math.pi/180
    root_chord = 2.0
    tip_chord = 1.0
    wing_avg_chord = 1.5
    
    m_LE = math.tan(le_sweep)
    AR = 5
    lamda = 0.5
    delta_m = 4*(lamda -1)/(AR*(lamda + 1))

    num_wings = 1
    num_wing_parts = [2]
    origin = [Vec3([0.0,0.0,0.0])]
    freestream_vel = 20.0

    wing_y = generate_spanwise_control_points(span_vortex_nodes)

    chord_distribution = [(1.0 - (1.0-lamda)*2.0*y)*root_chord for y in wing_y]
    camber = 0.0
    wing_twist = np.zeros(span_vortex_nodes)
    wing_sweep = np.zeros(span_vortex_nodes) # change it to a formula based on wing parameters

    num_rotors = 0
    num_blades = [0]
    elements = 48

    aircraft = Aircraft(num_rotors, num_wings)
    
    def build_wing_part(wp_idx):
        
        wing_part = WingPartGeometry(
            span_elements = span_vortex_nodes,
            chordwise_nodes = chord_vortex_nodes,
            wing_root_origin = origin[0],
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
            origin = origin[w_idx],
            wing_span = wing_span
        )

        wing.wing_parts = [build_wing_part(wp_idx) for wp_idx in range(num_wing_parts[w_idx])]
        
        return wing
    
    aircraft.wings= [build_wing(w_idx) for w_idx in range(num_wings)]
    set_wing_ctrl_pt_geometry(aircraft.wings[0], span_vortex_nodes, chord_vortex_nodes, camber)
    
    ac_state = AircraftState(num_rotors, num_blades, elements, num_wings, num_wing_parts, span_vortex_nodes, chord_vortex_nodes, aircraft)
    
    print("Freestream vel: ", freestream_vel, "m/s")
    
    ac_input_state = AircraftInputState(num_rotors, num_blades, num_wings)
    
    wing_lift_surface = WingLiftSurf(num_wing_parts[0])
    for wp_idx in range(num_wing_parts[0]):
        wing_part_lift_surf = WingPartLiftingSurf(span_vortex_nodes, chord_vortex_nodes)
        wing_lift_surface.wing_part_lift_surf[wp_idx] = wing_part_lift_surf
        
    set_wing_vortex_geometry(wing_lift_surface, aircraft.wings[0], span_vortex_nodes, chord_vortex_nodes)
        
    ac_input_state.wing_inputs[0].angle_of_attack = alpha
    ac_input_state.wing_inputs[0].freestream_velocity = freestream_vel
    
    # calculations 
    
    for wp_idx, wing_part_state in enumerate(ac_state.wing_states[0].wing_part_states):
        for sp_chunk in range(chord_vortex_nodes*span_vortex_nodes/8):
            for c1 in range(chunk_size()):
                wing_part_state.ctrl_chunks.ctrl_pt_aoa[c1] = alpha
                wing_part_state.ctrl_chunks.ctrl_pt_u_t[c1] = freestream_vel
    

                
    
    