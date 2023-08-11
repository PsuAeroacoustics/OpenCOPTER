
import sys
import os

sys.path.append(f'{os.path.dirname(os.path.realpath(__file__))}../../dependencies/OpenCOPTER')
sys.path.append(f'{os.path.dirname(os.path.realpath(__file__))}../../dependencies/wopwopd')

from libopencopter import *
from libwopwopd import *
import wopwop_input_files_generator
import simulate_aircraft
import numpy as np
import math
import scipy.io
from scipy.interpolate import interp1d
from scipy.misc import derivative
import numpy as np

from os import path, makedirs

def build_blade_from_json(blade_object, requested_elements, geom_directory):
    theta_tw_1 = blade_object['theta_tw']
    airfoil_descs = blade_object['airfoils']
    r_c = 0
    if "r_c" in blade_object:
        r_c = blade_object["r_c"]

    r = generate_radius_points(requested_elements, r_c)
    elements = len(r)

    non_dim_length = 1.0 - r_c
    
    airfoils = []
    extents = []
    for airfoil_desc in airfoil_descs:
        type = airfoil_desc['type']
        extent = [0, 0]
        af_extent = airfoil_desc['extent']
        airfoil = None
        if type == 'aerodas':
            airfoil = create_aerodas_from_xfoil_polar(f"{geom_directory}/{airfoil_desc['xfoil_polar']}", airfoil_desc['thickness'])
        elif type == 'thinaf':
            airfoil = ThinAirfoil(0.0401)
        else:
            print(f"Unsupported airfoil type: {type}. Defaulting to ThinAirfoil theory")
            airfoil = ThinAirfoil(0.0401)

        extent[0] = int(round((1.0 - (1.0/math.pi)*math.acos((2.0*(af_extent[0] - r_c)/non_dim_length) - 1.0))*elements))
        extent[1] = int(round((1.0 - (1.0/math.pi)*math.acos((2.0*(af_extent[1] - r_c)/non_dim_length) - 1.0))*elements)) - 1

        extents.append(extent)
        airfoils.append(airfoil)

    blade_airfoil = BladeAirfoil(airfoils, extents)

    x = np.zeros(elements)

    linear_x = None
    linear_r = None
    linear_c = None
    if 'x' in blade_object:
        linear_x = blade_object['x']
        linear_r = np.linspace(r_c, 1.0, len(linear_x))
        f_x = interp1d(linear_r, linear_x)
        x = f_x(r)

    c = None
    if "AR" in blade_object:
        AR = blade_object["AR"]
        c = ((1.0 - r_c)/AR)*np.ones(elements)
        if linear_r is not None:
            linear_c = ((1.0 - r_c)/AR)*np.ones(len(linear_r))
        else:
            linear_r = np.linspace(r_c, 1.0, elements)
            linear_x = np.zeros(len(linear_r))
            linear_c = ((1.0 - r_c)/AR)*np.ones(len(linear_r))

    else:
        linear_c = blade_object['c']
        if linear_r is None:
            linear_x = np.zeros(len(linear_c))
            linear_r = np.linspace(r_c, 1.0, len(linear_c))

        f_c = interp1d(linear_r, linear_c)
        c = f_c(r)

    x_over_c = np.asarray(linear_x)/np.asarray(linear_c)
    
    f_x_over_c = interp1d(linear_r, x_over_c, bounds_error=False, fill_value='extrapolate')

    linear_x_over_c_p = [derivative(f_x_over_c, _r, 1.0e-12) for _r in linear_r[3:-3]]

    f_x_over_c_p = interp1d(linear_r[3:-3], linear_x_over_c_p, bounds_error=False, fill_value='extrapolate')

    xp = f_x_over_c_p(r)

    twist = np.asarray([(_r - 0.75)*theta_tw_1*(math.pi/180.0) for _r in generate_radius_points(requested_elements, 0.0)])

    # Build the geom of the blades
    blade = BladeGeometry(
        num_elements = elements,
        azimuth_offset = 0,
        average_chord = 0,
        airfoil = blade_airfoil
    )
    
    sweep = sweep_from_quarter_chord(r, x)

    # Convert from linear array format to OpenCOPTER chunk format
    set_r(blade, r)
    set_xi(blade, x)
    set_xi_p(blade, xp)
    set_twist(blade, twist)
    set_chord(blade, c)
    set_sweep(blade, sweep)

    return blade, r_c

def compute_aero(log_file, args, V_inf, aoa, output_base, do_compute, geometry, flight_condition, computational_parameters, results, geom_directory):

    num_rotors = len(geometry["rotors"])

    omegas =  np.asarray([flight_condition["omegas"][r_idx] for r_idx in range(num_rotors)])
    d_psi = computational_parameters["d_psi"]*math.pi/180.0
    dt = d_psi/np.max(np.abs(omegas))

    requested_elements = computational_parameters["spanwise_elements"]

    

    R = np.asarray([geometry["rotors"][r_idx]["radius"] for r_idx in range(len(geometry["rotors"]))])

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

    sos = flight_condition["sos"]
    
    shed_history_angle = computational_parameters["shed_history_angle"]
    shed_history = round(shed_history_angle/(d_psi*180.0/math.pi))

    blade_dict = {}
    if "blades" in geometry:
        for blade_obj in geometry['blades']:
            blade_dict[blade_obj["name"]] = blade_obj

	# Create the spanwise distributions for our blades
	# We will later apply these to murosims internal
	# data layout. r and c arrays are normalized by the
	# rotor radius
    r = generate_radius_points(requested_elements)
    elements = len(r)
    log_file.write(f"requested_elements: {requested_elements}, actual elements:  {elements}\n")

    atmo = Atmosphere(density = density, dynamic_viscosity = dynamic_viscosity)
	
	# Create our outer rotorcraft_system geometry container
    rotorcraft_system = Aircraft(num_rotors)
	
	# The origins of our rotor
    origins = [Vec3(geometry['rotors'][r_idx]["origin"]) for r_idx in range(len(geometry['rotors']))]

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
    #wake_history_revs = np.round(wake_dists/(mus*2.0*math.pi)) # this is an array

    wake_history_revs = np.zeros(wake_dists.size)

    for idx, mu in enumerate(mus):
        if mu != 0:
            wake_history_revs[idx] = round(wake_dists[idx]/(mu*2.0*math.pi)) # this is an array
        else:
            wake_history_revs[idx] = round(wake_dists[idx]/(0.1*2.0*math.pi)) # this is an array

    log_file.write(f'wake_history_revs: {wake_history_revs}\n')

    wake_history_length = np.round(2*math.pi/d_psi*wake_history_revs).astype(int)
    wake_history_length = wake_history_length.tolist()

    naca0012_xsection = naca0012()

    wopwop_data_path = f'{output_base}/acoustics/data'

    if not path.isdir(wopwop_data_path):
        makedirs(wopwop_data_path, exist_ok=True)
                
    log_file.write(f'wake_history_length: {wake_history_length}\n') #this must be a list of intigers

    def build_blade(r_idx, b_idx, num_blades):

        blade = None
        r_c = 0
        if "number_of_blades" in geometry['rotors'][r_idx]:
            if isinstance(geometry['rotors'][r_idx]['blade'], str):
                blade, r_c = build_blade_from_json(blade_dict[geometry['rotors'][r_idx]['blade']], requested_elements, geom_directory)
            else:
                blade, r_c = build_blade_from_json(geometry['rotors'][r_idx]['blade'], requested_elements, geom_directory)
        else:
            if isinstance(geometry['rotors'][r_idx]['blades'], List[str]):
                blade, r_c = build_blade_from_json(blade_dict[geometry['rotors'][r_idx]['blade'][b_idx]], requested_elements, geom_directory)
            else:
                blade, r_c = build_blade_from_json(geometry['rotors'][r_idx]['blades'][b_idx], requested_elements, geom_directory)

        d_azimuth = 2.0*math.pi/num_blades

        blade.azimuth_offset = b_idx*d_azimuth

        c = get_chord(blade)
        blade.average_chord = R[r_idx]*np.sum(c)/len(c)
        blade.blade_length = R[r_idx]*(1.0 - r_c)

        wopwop_input_files_generator.write_wopwop_geometry(naca0012_xsection, r, get_twist(blade), R[r_idx], 10, wopwop_data_path, r_idx, b_idx)

        return blade

    def build_rotor(r_idx):
        num_blades = 0
        if "number_of_blades" in geometry['rotors'][r_idx]:
            num_blades = geometry['rotors'][r_idx]['number_of_blades']
        else:
            num_blades = len(geometry['rotors'][r_idx]['blades'])

        rotor = RotorGeometry(
			num_blades = num_blades,
			origin = origins[r_idx],
			radius = R[r_idx],
			solidity = 0 # we set this later because we don't have the average blade chord
		)

        tmp_blades = [build_blade(r_idx, b_idx, num_blades) for b_idx in range(num_blades)]

        rotor.blades = tmp_blades
        rotor.solidity = num_blades*rotor.blades[0].average_chord/(math.pi*rotor.radius)
        return rotor

    
    rotorcraft_system.rotors = [build_rotor(r_idx) for r_idx in range(num_rotors)]
    
    num_blades = [rotor.blades.length() for rotor in rotorcraft_system.rotors]

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
        if "aoa_rotors" in flight_condition:
            rotorcraft_input_state.rotor_inputs[r_idx].angle_of_attack = flight_condition["aoa_rotors"][r_idx]*(math.pi/180.0)
        else:
            rotorcraft_input_state.rotor_inputs[r_idx].angle_of_attack = aoa
        rotorcraft_input_state.rotor_inputs[r_idx].angular_velocity = omegas[r_idx]
        rotorcraft_input_state.rotor_inputs[r_idx].angular_accel = 0
        rotorcraft_input_state.rotor_inputs[r_idx].azimuth = azimuths[r_idx]
        rotorcraft_input_state.rotor_inputs[r_idx].freestream_velocity = V_inf
        for b_idx in range(num_blades[r_idx]):
            rotorcraft_input_state.rotor_inputs[r_idx].r_0[b_idx] = 0.1*rotorcraft_system.rotors[r_idx].blades[b_idx].average_chord/R[r_idx]
            rotorcraft_input_state.rotor_inputs[r_idx].blade_flapping[b_idx] = 0
            rotorcraft_input_state.rotor_inputs[r_idx].blade_flapping_rate[b_idx] = 0
            rotorcraft_input_state.rotor_inputs[r_idx].blade_pitches[b_idx] = collectives[r_idx]

    rotorcraft_inflows = [HuangPeters(4, 2, rotorcraft_system.rotors[r_idx], dt) if num_blades[r_idx] != 2 else HuangPeters(2, 2, rotorcraft_system.rotors[r_idx], dt) for r_idx in range(num_rotors)]
    #rotorcraft_inflows = [HuangPeters(2, 2, rotorcraft_system.rotors[r_idx], dt) for r_idx in range(num_rotors)]
    #rotorcraft_inflows = [HuangPeters(6, 4, rotorcraft_system.rotors[r_idx], dt) for r_idx in range(num_rotors)]
    #rotorcraft_inflows = [HuangPeters(8, 6, rotorcraft_system.rotors[r_idx], dt) for r_idx in range(num_rotors)]

	# Setup the wake history. We need at minimum 2 timesteps worth of history for the update.
	# Increasing the history increases computation time with the current implementation
    log_file.write(f'wake_history_length: {wake_history_length}\n')
    rotor_wake_history = WakeHistory(num_rotors, num_blades, wake_history_length, 2, elements, shed_history)

    sim_revs = wake_history_revs.max() + args.r
    log_file.write(f"sim_revs = {sim_revs}\n")

    wopwop_rotor_indx = 0
    
    rotorcraft_thrusts, rotorcraft_namelists = simulate_aircraft.simulate_aircraft(log_file, rotorcraft_system, rotorcraft_state, rotorcraft_input_state, rotorcraft_inflows, rotor_wake_history, atmo, omegas, V_inf, r, 1, elements, wopwop_rotor_indx, args.ws, f'{output_base}/vtu', f'{output_base}/acoustics', do_compute, flight_condition, computational_parameters["convergence_criteria"])

    if do_compute:
        results_dictionary = {}
        for r_idx in range(num_rotors):
            actual_wake_history = wake_history_length[r_idx] if wake_history_length[r_idx]%chunk_size() == 0 else wake_history_length[r_idx] + (chunk_size() - wake_history_length[r_idx]%chunk_size())
            wake_trajectories = np.zeros((num_blades[r_idx], 3, actual_wake_history))

            for b_idx in range(num_blades[r_idx]):
                wake_trajectories[b_idx, 0, :] = get_wake_x_component(rotor_wake_history.history[0].rotor_wakes[r_idx].tip_vortices[b_idx])
                wake_trajectories[b_idx, 1, :] = get_wake_y_component(rotor_wake_history.history[0].rotor_wakes[r_idx].tip_vortices[b_idx])
                wake_trajectories[b_idx, 2, :] = get_wake_z_component(rotor_wake_history.history[0].rotor_wakes[r_idx].tip_vortices[b_idx])

            results_dictionary[f'wake_{r_idx}_trajectory'] = wake_trajectories

        results_dictionary["rotor_c_t"] = rotorcraft_thrusts
        results_dictionary["rotor_collectives"] = [rotorcraft_input_state.rotor_inputs[r_idx].blade_pitches[0] for r_idx in range(num_rotors)]
        results_dictionary["rotor_chis"] = [rotorcraft_inflows[r_idx].wake_skew() for r_idx in range(num_rotors)]

        scipy.io.savemat(f"{output_base}/results.mat", results_dictionary)

        for slice_idx, inflow_slice in enumerate(results["inflow_slices"]):
            res_x = inflow_slice["resolution"][0]
            res_y = inflow_slice["resolution"][1]
            res_z = inflow_slice["resolution"][2]

            deltas = Vec3([inflow_slice["slice_size"][0]/res_x, inflow_slice["slice_size"][1]/res_y, inflow_slice["slice_size"][2]/res_z])
            start = Vec3(inflow_slice["slice_start"])

            write_inflow_vtu(f"{output_base}/inflow_model_slice_{slice_idx}.vtu", rotorcraft_inflows, deltas, start, res_x, res_y, res_z, origins, rotorcraft_input_state.rotor_inputs[0].angle_of_attack, omegas[0])

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
