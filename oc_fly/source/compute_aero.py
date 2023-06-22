
import sys
import os

sys.path.append(f'{os.path.dirname(os.path.realpath(__file__))}../../dependencies/OpenCOPTER')
sys.path.append(f'{os.path.dirname(os.path.realpath(__file__))}../../dependencies/wopwopd')

from libopencopter import *
from libwopwopd import *
import simulate_aircraft
import numpy as np
import math
import scipy.io
from os import path, makedirs

def compute_aero(log_file, args, V_inf, aoa, output_base, do_compute, geometry, flight_condition, computational_parameters):

    aoa_rotors = flight_condition["aoa_rotors"]   # angle of attack for each rotor rather than for entire aircraft

    omegas =  np.asarray(flight_condition["omegas"])
    d_psi = computational_parameters["d_psi"]*math.pi/180.0
    dt = d_psi/np.max(np.abs(omegas))

    requested_elements = computational_parameters["spanwise_elements"]

    num_rotors = geometry["number_of_rotors"]
    num_blades = geometry["number_of_blades"]
    R = np.asarray(geometry["radii"])

    mus = V_inf/(np.abs(omegas)*R)

    log_file.write(f'mus: {mus}\n')
    log_file.write(f'R: {R}\n')
    log_file.write(f'omegas: {omegas}\n')

    collectives = None
    if "collective" in flight_condition:
        collectives = np.asarray(flight_condition["collective"])*np.pi/180
    else:
        collectives = 2.0*np.ones(num_rotors)*np.pi/180

    density = flight_condition["density"]
    dynamic_viscosity = 18.03e-6

    AR = geometry['AR']
    sos = flight_condition["sos"]
    
    #V_inf = param["flight conditions"][0]["V inf"]
    theta_tw_1 = geometry["theta_tw_1"]*math.pi/180.0

    shed_history_angle = computational_parameters["shed_history_angle"]
    shed_history = round(shed_history_angle/(d_psi*180.0/math.pi))

    d_azimuth = 2.0*math.pi/num_blades

	# Create the spanwise distributions for our blades
	# We will later apply these to murosims internal
	# data layout. r and c arrays are normalized by the
	# rotor radius
    r = generate_radius_points(requested_elements)
    elements = len(r)
    log_file.write(f"requested_elements: {requested_elements}, actual elements:  {elements}\n")
    c = (1.0/AR)*np.ones(elements)
    alpha_0 = -2.5*(math.pi/180.0)*np.ones(elements)
    twist = [(_r - 0.75)*theta_tw_1 for _r in r]
    C_l_alpha = 2.0*math.pi*np.ones(elements)
    sweep = np.zeros(elements)

    atmo = Atmosphere(density = density, dynamic_viscosity = dynamic_viscosity)
	
	# Create our outer rotorcraft_system geometry container
    rotorcraft_system = Aircraft(num_rotors)
	
	# The origins of our rotor
    origins = geometry["origin"]
    origins = [Vec3(origin) for origin in origins]

    wake_dists = np.zeros(num_rotors)

    for o_idx, origin in enumerate(origins):
        for other_origin in origins:
            delta_x = origin[0] - other_origin[0]

            # plus 1R to finish going past the other rotor, and then 6 more radii after that
            wake_dists[o_idx] = max(wake_dists[o_idx], delta_x + 1 + 6)

        if wake_dists[o_idx] <= 0:
            # if all others are upstream, just trail 6 radii behind
            wake_dists[o_idx] = 6

    log_file.write(f'wake_dists: {wake_dists}\n')
    wake_history_revs = np.round(wake_dists/(mus*2.0*math.pi)) # this is an array

    log_file.write(f'wake_history_revs: {wake_history_revs}\n')
    #log_file.write(type(wake_history_revs))

    wake_history_length = np.round(2*math.pi/d_psi*wake_history_revs).astype(int)
    wake_history_length = wake_history_length.tolist()
    #wake_history_length[0] = 2

    log_file.write(f'wake_history_length: {wake_history_length}\n') #this must be a list of intigers

    def build_blade(r_idx, b_idx):
		# Build the geom of the blades
        blade = BladeGeometry(
            num_elements = elements,
            azimuth_offset = b_idx*d_azimuth,
            average_chord = R[r_idx]*np.sum(c)/len(c)
        )
        
        # Convert from linear array format to murosim chunk format
        set_r(blade, r)
        set_twist(blade, twist)
        set_alpha_0(blade, alpha_0)
        set_C_l_alpha(blade, C_l_alpha)
        set_chord(blade, c)
        set_sweep(blade, sweep)

        return blade

    def build_rotor(r_idx):
        rotor = RotorGeometry(
			num_blades = num_blades,
			origin = origins[r_idx],
			radius = R[r_idx],
			solidity = 0 # we set this later because we don't have the average blade chord
		)

        rotor.blades = [build_blade(r_idx, b_idx) for b_idx in range(num_blades)]
        rotor.solidity = num_blades*rotor.blades[0].average_chord/(math.pi*rotor.radius)
        return rotor

    rotorcraft_system.rotors = [build_rotor(r_idx) for r_idx in range(num_rotors)]

	# AircraftState is the top level container for holding the current
	# aerodynamic state of the tandem_system. It breaks down into rotors
	# the the individual blades. There is a series of functions provided
	# to turn internal state data into a linear array.
    rotorcraft_state = None
    if do_compute:
        rotorcraft_state = AircraftState(num_rotors, num_blades, elements, rotorcraft_system)

    log_file.write(f"Freestream vel: {V_inf} m/s\n")

    azimuths = [0, 0*(math.pi/180.0)]
	# Create and setup the input state. This would be the
	# sort of input a dynamics simulator might feed into
	# the aero model.
    rotorcraft_input_state = AircraftInputState(num_rotors, num_blades)

    for r_idx in range(num_rotors):
        if aoa_rotors:
            rotorcraft_input_state.rotor_inputs[r_idx].angle_of_attack = aoa_rotors[r_idx]*(math.pi/180.0)
        else:
            rotorcraft_input_state.rotor_inputs[r_idx].angle_of_attack = aoa
        rotorcraft_input_state.rotor_inputs[r_idx].angular_velocity = omegas[r_idx]
        rotorcraft_input_state.rotor_inputs[r_idx].angular_accel = 0
        rotorcraft_input_state.rotor_inputs[r_idx].azimuth = azimuths[r_idx]
        rotorcraft_input_state.rotor_inputs[r_idx].freestream_velocity = V_inf
        for b_idx in range(num_blades):
            rotorcraft_input_state.rotor_inputs[r_idx].r_0[b_idx] = 0.1*rotorcraft_system.rotors[r_idx].blades[b_idx].average_chord
            rotorcraft_input_state.rotor_inputs[r_idx].blade_flapping[b_idx] = 0
            rotorcraft_input_state.rotor_inputs[r_idx].blade_flapping_rate[b_idx] = 0
            rotorcraft_input_state.rotor_inputs[r_idx].blade_pitches[b_idx] = collectives[r_idx]

    #rotorcraft_inflows = [HuangPeters(4, 2, rotorcraft_system.rotors[r_idx], dt) if r_idx > 0 else Decayed() for r_idx in range(num_rotors)]
    rotorcraft_inflows = [HuangPeters(4, 2, rotorcraft_system.rotors[r_idx], dt) for r_idx in range(num_rotors)]

	# Setup the wake history. We need at minimum 2 timesteps worth of history for the update.
	# Increasing the history increases computation time with the current implementation
    log_file.write(f'wake_history_length: {wake_history_length}\n')
    rotor_wake_history = WakeHistory(num_rotors, num_blades,wake_history_length, 2, elements, shed_history)

    sim_revs = wake_history_revs.max() + args.r
    log_file.write(f"sim_revs = {sim_revs}\n")

    wopwop_rotor_indx = 0
    # if (num_rotors == 1):
    #     wopwop_rotor_indx = 0
    # else:
    #     wopwop_rotor_indx = computational_parameters["wopwop_rotor_indx"]

    #rotorcraft_thrusts, rotorcraft_namelists, x_inflow, inflow_induced_velocities, wake_induced_velocities = simulate_aircraft.simulate_aircraft(rotorcraft_system, rotorcraft_state, rotorcraft_input_state, rotorcraft_inflows, rotor_wake_history, atmo, omegas, V_inf, r, twist, AR, 1, elements, wopwop_rotor_indx, args.ws, f'{output_base}/vtu', f'{output_base}/acoustics', do_compute, flight_condition, computational_parameters["convergence_criteria"])
    rotorcraft_thrusts, rotorcraft_namelists = simulate_aircraft.simulate_aircraft(log_file, rotorcraft_system, rotorcraft_state, rotorcraft_input_state, rotorcraft_inflows, rotor_wake_history, atmo, omegas, V_inf, r, twist, AR, 1, elements, wopwop_rotor_indx, args.ws, f'{output_base}/vtu', f'{output_base}/acoustics', do_compute, flight_condition, computational_parameters["convergence_criteria"])

    if do_compute:
        results_dictionary = {}
        for r_idx in range(num_rotors):
            actual_wake_history = wake_history_length[r_idx] if wake_history_length[r_idx]%chunk_size() == 0 else wake_history_length[r_idx] + (chunk_size() - wake_history_length[r_idx]%chunk_size())
            wake_trajectories = np.zeros((num_blades, 3, actual_wake_history))

            for b_idx in range(num_blades):
                wake_trajectories[b_idx, 0, :] = get_wake_x_component(rotor_wake_history.history[0].rotor_wakes[r_idx].tip_vortices[b_idx])
                wake_trajectories[b_idx, 1, :] = get_wake_y_component(rotor_wake_history.history[0].rotor_wakes[r_idx].tip_vortices[b_idx])
                wake_trajectories[b_idx, 2, :] = get_wake_z_component(rotor_wake_history.history[0].rotor_wakes[r_idx].tip_vortices[b_idx])

            results_dictionary[f'wake_{r_idx}_trajectory'] = wake_trajectories

        results_dictionary["rotor_c_t"] = rotorcraft_thrusts
        results_dictionary["rotor_collectives"] = [rotorcraft_input_state.rotor_inputs[r_idx].blade_pitches[0] for r_idx in range(num_rotors)]
        results_dictionary["rotor_chis"] = [rotorcraft_inflows[r_idx].wake_skew() for r_idx in range(num_rotors)]
    
        #results_dictionary["x_inflow"] = x_inflow
        #results_dictionary["wake_induced_velocities"] = wake_induced_velocities
        #results_dictionary["inflow_induced_velocities"] = inflow_induced_velocities
        
        scipy.io.savemat(f"{output_base}/results.mat", results_dictionary)
    
    cases = []
    for r_idx, namelist in enumerate(rotorcraft_namelists[0:num_rotors]):
        wopwop_case_path = f'{output_base}/acoustics/rotor_{r_idx}/'

        if not path.isdir(wopwop_case_path):
            makedirs(wopwop_case_path, exist_ok=True)

        single_rotor_case = Casename()
        single_rotor_case.caseNameFile = "case.nam"
        single_rotor_case.globalFolderName = wopwop_case_path
        single_rotor_case.namelist = namelist

        cases.append(single_rotor_case)


    wopwop_case_path = f'{output_base}/acoustics/full_system/'

    if not path.isdir(wopwop_case_path):
        makedirs(wopwop_case_path, exist_ok=True)

    full_system_case = Casename()
    full_system_case.caseNameFile = "case.nam"
    full_system_case.globalFolderName = wopwop_case_path
    full_system_case.namelist = rotorcraft_namelists[-1]

    cases.append(full_system_case)

    return cases
