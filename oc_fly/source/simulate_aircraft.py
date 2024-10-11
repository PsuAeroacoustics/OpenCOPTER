from ast import operator
import sys
import os

sys.path.append(f'{os.path.dirname(os.path.realpath(__file__))}/dependencies/OpenCOPTER')
sys.path.append(f'{os.path.dirname(os.path.realpath(__file__))}/dependencies/wopwopd')

from libopencopter import *
from libwopwopd import *
import wopwop_input_files_generator
import numpy as np
from scipy.integrate import simpson
import math
import time

from simulated_vehicle import SimulatedVehicle

from os import path, makedirs

def flapping_at_azimuth(a: list[float], b: list[float], w: float, azimuth: float):
	h = 0
	h_star = 0

	for idx in range(len(a)):
		cos = math.cos(w*float(idx)*azimuth)
		sin = math.sin(w*float(idx)*azimuth)

		h = h + a[idx]*cos + b[idx - 1]*sin

		if idx > 0:
			h_star = h_star + (-w*float(idx)*a[idx]*sin) + w*float(idx)*b[idx - 1]*cos

	return (h, h_star)

def elastic_twist_at_azimuth(a: list[float], b: list[float], azimuth: float):
	h = 0
	#h_star = 0

	for idx in range(len(a)):
		cos = math.cos(float(idx)*azimuth)
		sin = math.sin(float(idx)*azimuth)

		h = h + a[idx]*cos + b[idx - 1]*sin

	return h

def simulate_aircraft(log_file, vehicle: SimulatedVehicle, atmo, elements, write_wake, vtu_output_path, wopwop_output_path, do_compute, flight_condition, computational_parameters, observer, acoustics, wake_lengths, results, wopwop_motion):

	if not path.isdir(wopwop_output_path):
		makedirs(wopwop_output_path, exist_ok=True)

	if not path.isdir(vtu_output_path):
		makedirs(vtu_output_path, exist_ok=True)

	num_rotors = vehicle.input_state.rotor_inputs.length()

	omegas = np.asarray([vehicle.input_state.rotor_inputs[rotor_idx].angular_velocity for rotor_idx in range(num_rotors)])

	log_file.write(f'num_rotors: {num_rotors}\n')
	num_blades = [vehicle.aircraft.rotors[rotor_idx].blades.length() for rotor_idx in range(vehicle.aircraft.rotors.length())]

	d_psi = computational_parameters['d_psi']

	dt = d_psi*(math.pi/180.0)/np.max(np.abs(omegas))
	iter_per_rev = 360/d_psi

	vtk_rotors = [build_base_vtu_rotor(vehicle.aircraft.rotors[rotor_idx]) for rotor_idx in range(num_rotors)]
	vtk_wake = build_base_vtu_wake(vehicle.wake_history.history[0])

	#C_T_len = int(round(2.0*math.pi/(dt*max(abs(omegas)))))
	C_T_len = np.round(2.0*math.pi/(dt*np.abs(omegas))).astype(dtype=np.int64)

	num_rotors = vehicle.input_state.rotor_inputs.length()
	average_C_T_arrays = [np.zeros(C_T_len[r_idx]) for r_idx in range(num_rotors)]
	average_C_Mx_arrays = [np.zeros(C_T_len[r_idx]) for r_idx in range(num_rotors)]
	average_C_My_arrays = [np.zeros(C_T_len[r_idx]) for r_idx in range(num_rotors)]
	average_C_Q_arrays = [np.zeros(C_T_len[r_idx]) for r_idx in range(num_rotors)]
	average_theta_arrays = [np.zeros(C_T_len[r_idx]) for r_idx in range(num_rotors)]
	average_theta_1c_arrays = [np.zeros(C_T_len[r_idx]) for r_idx in range(num_rotors)]
	average_theta_1s_arrays = [np.zeros(C_T_len[r_idx]) for r_idx in range(num_rotors)]

	average_C_Ts = np.zeros(num_rotors)
	average_C_Mxs = np.zeros(num_rotors)
	average_C_Mys = np.zeros(num_rotors)
	average_C_Qs = np.zeros(num_rotors)
	average_Qs = np.zeros(num_rotors)
	average_Ts = np.zeros(num_rotors)
	average_Mxs = np.zeros(num_rotors)
	average_Mys = np.zeros(num_rotors)
	average_theta = np.zeros(num_rotors)
	average_theta_1c = np.zeros(num_rotors)
	average_theta_1s = np.zeros(num_rotors)

	average_chi_arrays = [np.zeros(C_T_len[r_idx]) for r_idx in range(num_rotors)]
	average_chis = np.zeros(num_rotors)

	loading_data = ZoneVectorData(elements)

	loading_data.set_x_loading_array(np.zeros(elements, dtype=np.single))
	loading_data.set_y_loading_array(np.zeros(elements, dtype=np.single))

	naca0012_xsection = naca0012()

	wopwop_data_path = f'{wopwop_output_path}/data'

	loading_files = []

	post_conv_revolutions = 2
	if "post_conv_revolutions" in computational_parameters:
		post_conv_revolutions = computational_parameters["post_conv_revolutions"]

	log_file.write(f"post_conv_revolutions: {post_conv_revolutions}\n")
	if do_compute:
		loading_files = [[wopwop_input_files_generator.build_wopwop_loading(rotor, blade, int(round(post_conv_revolutions*iter_per_rev)), naca0012_xsection, wopwop_data_path) for blade in rotor.blades] for rotor in vehicle.aircraft.rotors]
	
	for rotor in vehicle.aircraft.rotors:
		for blade in rotor.blades:
			wopwop_input_files_generator.write_wopwop_geometry(naca0012_xsection, wopwop_data_path, rotor, blade, acoustics["thickness_noise_flag"])

	log_file.write("Performing acoustic resolution simulation\n")

	trim = False
	moment_trim = False

	trim_algo = 'he'

	c_t_bars = None
	c_mx_bars = None
	c_my_bars = None
	if "c_t" in flight_condition:
		trim = True
		c_t_bars = flight_condition["c_t"]

	elif "T" in flight_condition:
		trim = True
		c_t_bars =  [flight_condition["T"][rotor_idx]/(atmo.density*math.pi*vehicle.aircraft.rotors[rotor_idx].radius**4.0*abs(omegas[rotor_idx])**2.0) for rotor_idx in range(num_rotors)]

	if "c_mx" in flight_condition and "c_my" in flight_condition:
		moment_trim = True
		c_mx_bars = [flight_condition["c_mx"][rotor_idx]*atmo.density*math.pi*vehicle.aircraft.rotors[rotor_idx].radius**5.0*abs(omegas[rotor_idx])**2.0 for rotor_idx in range(num_rotors)]
		c_my_bars = [flight_condition["c_my"][rotor_idx]*atmo.density*math.pi*vehicle.aircraft.rotors[rotor_idx].radius**5.0*abs(omegas[rotor_idx])**2.0 for rotor_idx in range(num_rotors)]

	elif "Mx" in flight_condition and "My" in flight_condition:
		moment_trim = True
		c_mx_bars = flight_condition["Mx"]
		c_my_bars = flight_condition["My"]

		log_file.write(f'c_mx_bars: {c_mx_bars}, c_my_bars: {c_my_bars}\n')

	if "trim_algo" in computational_parameters:
		trim_algo = computational_parameters["trim_algo"]

	Ks = 0.4*np.ones(num_rotors)
	taus = 4*np.ones(num_rotors)
	moment_Ks = np.ones((num_rotors, 2))
	moment_taus = np.ones((num_rotors, 2))

	# HART II aerodas constants
	moment_Ks[0, 0] = 0.18
	moment_Ks[0, 1] = 0.18

	moment_taus[0, 0] = .5
	moment_taus[0, 1] = .5

	thetas = np.zeros((num_rotors, 2))
	last_thetas = np.zeros(num_rotors)
	moment_thetas = np.zeros((num_rotors, 2*2))

	theta_1s = np.zeros(num_rotors)
	theta_1c = np.zeros(num_rotors)

	for rotor_idx in range(num_rotors):
		last_thetas[rotor_idx] = 10000
		thetas[rotor_idx, 1] = vehicle.input_state.rotor_inputs[rotor_idx].blade_pitches[0]

	curr_c_ts = np.zeros(num_rotors)
	curr_c_mxs = np.zeros(num_rotors)
	curr_c_mys = np.zeros(num_rotors)

	last_wake_points = [None for _ in range(num_rotors)]
	last_C_T = [100 for _ in range(num_rotors)]
	last_C_Mx = [100 for _ in range(num_rotors)]
	last_C_My = [100 for _ in range(num_rotors)]

	wake_l2 = [1000 for _ in range(num_rotors)]
	sim_done = False
	converged = False
	converged_revolutions = 0
	iteration = 0
	acoustic_iteration = 0
	spanwise_element_iteration = [0 for _ in range(num_rotors)]
	aoa_update_iter = 0

	start_time = time.perf_counter_ns()

	theta_3 = 0
	psi_3 = 0

	temp_wake_array_x = np.zeros(wake_lengths[0])
	temp_wake_array_y = np.zeros(wake_lengths[0])
	temp_wake_array_z = np.zeros(wake_lengths[0])
	temp_wake_array_r_c = np.zeros(wake_lengths[0])

	track_wake_element = 'element_trajectories' in results
	track_span_element = 'spanwise_time_series' in results

	start_recording = False
	done_recording = False

	target_span_elements = []
	target_span_chunk_index = []
	target_span_element_index = []

	if track_span_element:
		target_span_elements = results['spanwise_time_series']

		target_span_chunk_index = [[0 for _ in range(len(target_span_elements))] for _ in range(num_rotors)]
		target_span_element_index = [[0 for _ in range(len(target_span_elements))] for _ in range(num_rotors)]

		for t_idx, target_span_element in enumerate(target_span_elements):
			closest_element_dist = math.inf
			for rotor_idx in range(num_rotors):
				for chunk_idx in range(vehicle.aircraft.rotors[rotor_idx].blades[0].chunks.len()):
					for sub_idx in range(chunk_size()):
						element_dist = abs(target_span_element - vehicle.aircraft.rotors[rotor_idx].blades[0].chunks[chunk_idx].r[sub_idx])
						if element_dist < closest_element_dist:
							target_span_chunk_index[rotor_idx][t_idx] = chunk_idx
							target_span_element_index[rotor_idx][t_idx] = sub_idx
							closest_element_dist = element_dist

	collective_pitch_array = np.zeros((num_rotors, int(round(post_conv_revolutions*iter_per_rev)) + 1))
	cos_pitch_array = np.zeros((num_rotors, max(num_blades), int(round(post_conv_revolutions*iter_per_rev)) + 1))
	sin_pitch_array = np.zeros((num_rotors, max(num_blades), int(round(post_conv_revolutions*iter_per_rev)) + 1))
	hhc_pitch_array = np.zeros((num_rotors, max(num_blades), int(round(post_conv_revolutions*iter_per_rev)) + 1))
	blade_flapping_array = np.zeros((num_rotors, max(num_blades), int(round(post_conv_revolutions*iter_per_rev)) + 1))
	blade_flapping_der_array = np.zeros((num_rotors, max(num_blades), int(round(post_conv_revolutions*iter_per_rev)) + 1))
	elastic_twist_array = np.zeros((num_rotors, max(num_blades), int(round(post_conv_revolutions*iter_per_rev)) + 1))
	blade_twist_array = np.zeros((num_rotors, max(num_blades), int(round(post_conv_revolutions*iter_per_rev)) + 1))
	blade_twist_azimuth = np.zeros((num_rotors, max(num_blades), int(round(post_conv_revolutions*iter_per_rev)) + 1))

	# Nitya: checking if these gets stored in the matlab file!

	wake_idx =  [[[[] for _ in range(int(round(post_conv_revolutions*iter_per_rev)) + 1)] for _ in range(max(num_blades))] for _ in range(max(num_blades))]
	#bladeSec_idx = [[[[] for _ in range(int(round(post_conv_revolutions*iter_per_rev)) + 1)] for _ in range(max(num_blades))] for _ in range(max(num_blades))]
	#wake_miss_dist = [[[[] for _ in range(int(round(post_conv_revolutions*iter_per_rev)) + 1)] for _ in range(max(num_blades))] for _ in range(max(num_blades))]
	blade_directionVec = [[[[] for _ in range(int(round(post_conv_revolutions*iter_per_rev)) + 1)] for _ in range(max(num_blades))] for _ in range(max(num_blades))]
	vortex_directionVec = [[[[] for _ in range(int(round(post_conv_revolutions*iter_per_rev)) + 1)] for _ in range(max(num_blades))] for _ in range(max(num_blades))]

	target_y_slices = []
	if track_wake_element:
		target_y_slices = results['element_trajectories']

	span_element_loading = np.zeros((num_rotors, len(target_span_elements), int(post_conv_revolutions*iter_per_rev)))
	span_element_aoa_eff = np.zeros((num_rotors, len(target_span_elements), int(post_conv_revolutions*iter_per_rev)))
	span_element_aoa = np.zeros((num_rotors, len(target_span_elements), int(post_conv_revolutions*iter_per_rev)))
	span_element_up = np.zeros((num_rotors, len(target_span_elements), int(post_conv_revolutions*iter_per_rev)))
	span_element_inflow_angle = np.zeros((num_rotors, len(target_span_elements), int(post_conv_revolutions*iter_per_rev)))
	span_element_theta = np.zeros((num_rotors, len(target_span_elements), int(post_conv_revolutions*iter_per_rev)))
	span_element_gamma = np.zeros((num_rotors, len(target_span_elements), int(post_conv_revolutions*iter_per_rev)))

	wake_element_index = np.zeros((num_rotors, len(target_y_slices)), dtype=int)
	wake_element_blade = np.zeros((num_rotors, len(target_y_slices)), dtype=int)
	wake_element_found = [[False for _ in range(len(target_y_slices))] for _ in range(num_rotors)]

	wake_element_trajectory = np.zeros((num_rotors, len(target_y_slices), 2, int(post_conv_revolutions*iter_per_rev) + 1))
	wake_element_core_size = np.zeros((num_rotors, len(target_y_slices), int(post_conv_revolutions*iter_per_rev) + 1))
	
	convergence_type = 'wake'
	if 'convergence_type' in computational_parameters:
		convergence_type = computational_parameters['convergence_type']
		if convergence_type == 'wake':
			convergence_orientation = 'z'
			if 'wake_convergence_orientation' in computational_parameters:
				convergence_orientation = computational_parameters['wake_convergence_orientation']

	if 'psi_3' in flight_condition:
		psi_3 = flight_condition['psi_3']*(math.pi/180.0)
	
	if 'theta_3' in flight_condition:
		theta_3 = flight_condition['theta_3']*(math.pi/180.0)

	blade_flapping = None
	elastic_twist = None

	# for rotor_idx in range(num_rotors):
	# 	write_rotor_vtu(f"{vtu_output_path}/rotor", 1000000000, rotor_idx, vtk_rotors[rotor_idx], vehicle.ac_state.rotor_states[rotor_idx], vehicle.input_state.rotor_inputs[rotor_idx], vehicle.aircraft.rotors[rotor_idx])
	
	convergence_rev_multiple = 1
	if "convergence_rev_multiple" in computational_parameters:
		convergence_rev_multiple = computational_parameters["convergence_rev_multiple"]

	max_l2 = 1000

	if do_compute:

		for m in flight_condition["motion"]:
			if 'blade_element_func' in m:
				if m['blade_element_func'] == 'flapping':
					blade_flapping = lambda az: flapping_at_azimuth(m["cos"], m["sin"], 1.0, az)
				elif m['blade_element_func'] == 'pitching':
					elastic_twist = lambda az: elastic_twist_at_azimuth(m["cos"], m["sin"], az)

		z_loading = np.zeros(elements, dtype=np.single)
		x_loading = np.zeros(elements, dtype=np.single)
		y_loading = np.zeros(elements, dtype=np.single)

		if convergence_type == 'wake':
			for rotor_idx in range(num_rotors):
				if convergence_orientation == 'x':
					last_wake_points[rotor_idx] = np.asarray(get_wake_x_component(vehicle.wake_history.history[0].rotor_wakes[rotor_idx].tip_vortices[0]))
				elif convergence_orientation == 'y':
					last_wake_points[rotor_idx] = np.asarray(get_wake_y_component(vehicle.wake_history.history[0].rotor_wakes[rotor_idx].tip_vortices[0]))
				elif convergence_orientation == 'z':
					last_wake_points[rotor_idx] = np.asarray(get_wake_z_component(vehicle.wake_history.history[0].rotor_wakes[rotor_idx].tip_vortices[0]))

		while not sim_done:

			
			if (iteration > 0) and (iteration % int(convergence_rev_multiple*iter_per_rev) == 0):
				max_l2 = 1000
				#log_file.write(f'checking convergence itr: {iteration}, convergence_rev_multiple*iter_per_rev: {convergence_rev_multiple*iter_per_rev}\n')
				for rotor_idx in range(num_rotors):
					if convergence_type == 'wake':
						wake_points = []
						if convergence_orientation == 'x':
							wake_points = np.asarray(get_wake_x_component(vehicle.wake_history.history[0].rotor_wakes[rotor_idx].tip_vortices[0]))
						elif convergence_orientation == 'y':
							wake_points = np.asarray(get_wake_y_component(vehicle.wake_history.history[0].rotor_wakes[rotor_idx].tip_vortices[0]))
						elif convergence_orientation == 'z':
							wake_points = np.asarray(get_wake_z_component(vehicle.wake_history.history[0].rotor_wakes[rotor_idx].tip_vortices[0]))

						wake_l2[rotor_idx] = np.sqrt(np.mean(np.power((wake_points - last_wake_points[rotor_idx]), 2.0)))

						last_wake_points[rotor_idx] = wake_points

						if "convergence_rotor" in computational_parameters:
							max_l2 = wake_l2[computational_parameters["convergence_rotor"]]
						else:
							max_l2 = np.abs(wake_l2).max()

					elif convergence_type == 'trim':
						C_T_delta = np.abs(average_C_T_arrays[rotor_idx][iteration % C_T_len[rotor_idx]] - last_C_T[rotor_idx])
						C_Mx_delta = np.abs(average_C_Mx_arrays[rotor_idx][iteration % C_T_len[rotor_idx]] - last_C_Mx[rotor_idx])
						C_My_delta = np.abs(average_C_My_arrays[rotor_idx][iteration % C_T_len[rotor_idx]] - last_C_My[rotor_idx])

						max_l2 = max([C_T_delta, C_Mx_delta, C_My_delta])

						last_C_T[rotor_idx] = average_C_T_arrays[rotor_idx][iteration % C_T_len[rotor_idx]]
						last_C_Mx[rotor_idx] = average_C_Mx_arrays[rotor_idx][iteration % C_T_len[rotor_idx]]
						last_C_My[rotor_idx] = average_C_My_arrays[rotor_idx][iteration % C_T_len[rotor_idx]]

				if (max_l2 is not None) and (max_l2 <= computational_parameters["convergence_criteria"]):
					if not converged:
						log_file.write(f"Simulation reached convergence criteria: {max_l2}\n")

					converged = True

			#if iteration % 1 == 0:
			if iteration % iter_per_rev == 0:
				now = time.perf_counter_ns()
				elapsed = now - start_time

				log_file.write(
					f'{elapsed/(1000**3):.5f}: rotor rev: {iteration/iter_per_rev:.3f},'
					+''.join([f' C_T{rotor_idx}: {average_C_Ts[rotor_idx]:.4f},' for rotor_idx in range(num_rotors)])
					+''.join([f' T{rotor_idx}: {average_Ts[rotor_idx]:.4f},' for rotor_idx in range(num_rotors)])
					+''.join([f' P{rotor_idx}: {average_Qs[rotor_idx]:.4f},' for rotor_idx in range(num_rotors)])
					+''.join([f' Mx_{rotor_idx}: {average_Mxs[rotor_idx]:.4f},' for rotor_idx in range(num_rotors)])
					+''.join([f' My_{rotor_idx}: {average_Mys[rotor_idx]:.4f},' for rotor_idx in range(num_rotors)])
					+''.join([f' θ{rotor_idx}: {average_theta[rotor_idx]*(180.0/math.pi):.4f},' for rotor_idx in range(num_rotors)])
					+''.join([f' θ1c_{rotor_idx}: {average_theta_1c[rotor_idx]*(180.0/math.pi):.4f},' for rotor_idx in range(num_rotors)])
					+''.join([f' θ1s_{rotor_idx}: {average_theta_1s[rotor_idx]*(180.0/math.pi):.4f},' for rotor_idx in range(num_rotors)])
					+''.join([f' χ{rotor_idx}: {average_chis[rotor_idx]:.2f},' for rotor_idx in range(num_rotors)])
					+f' combined_C_T: {average_C_Ts.sum():.5f}, max L_2: {max_l2}\n'
				)
				# log_file.write(f'theta: {[theta*(180.0/math.pi) for theta in get_theta(vehicle.ac_state.rotor_states[1].blade_states[0])]}\n')
				# log_file.write(f'aoa: {[aoa*(180.0/math.pi) for aoa in get_aoa(vehicle.ac_state.rotor_states[1].blade_states[0])]}\n')
				# log_file.write(f'aoa_eff: {[aoa*(180.0/math.pi) for aoa in get_aoa_eff(vehicle.ac_state.rotor_states[1].blade_states[0])]}\n')
				# log_file.write(f'inflow_angle: {[aoa*(180.0/math.pi) for aoa in get_inflow_angle(vehicle.ac_state.rotor_states[1].blade_states[0])]}\n')
				# log_file.write(f'u_p: {get_u_p(vehicle.ac_state.rotor_states[1].blade_states[0])}\n')
				# log_file.write(f'u_t: {get_u_t(vehicle.ac_state.rotor_states[1].blade_states[0])}\n')
				# log_file.write(f'gamma: {get_gamma(vehicle.ac_state.rotor_states[1].blade_states[0])}\n')
				# log_file.write(f'd_gamma: {get_d_gamma(vehicle.ac_state.rotor_states[1].blade_states[0])}\n')
				# log_file.write(f'dC_T: {get_dC_T(vehicle.ac_state.rotor_states[1].blade_states[0])}\n')

				start_time = now
				log_file.flush()

				if converged and not sim_done:
					if converged_revolutions >= post_conv_revolutions:
						sim_done = True

					converged_revolutions = converged_revolutions + 1

			for rotor_idx, rotor in enumerate(vehicle.aircraft.rotors):
				_ = basic_single_rotor_dynamics(vehicle.input_state.rotor_inputs[rotor_idx], dt)

				for motion_lambda in vehicle.motion_lambdas[rotor_idx]:
					motion_lambda(vehicle.input_state.rotor_inputs[rotor_idx].azimuth)
					# Nitya: The value of a is taken from here!! This will plug in the current azimuth value into all those functions defined in motion_lambda

					# Nitya: Where are we using these motion_lambdas?? 

			step(vehicle.ac_state, vehicle.aircraft, vehicle.input_state, vehicle.inflows, vehicle.wake_history, atmo, iteration, dt, converged)

			for rotor_idx, rotor_state in enumerate(vehicle.ac_state.rotor_states):
				average_C_T_arrays[rotor_idx][iteration % C_T_len[rotor_idx]] = rotor_state.C_T
				average_C_Q_arrays[rotor_idx][iteration % C_T_len[rotor_idx]] = rotor_state.C_Q
				average_C_Mx_arrays[rotor_idx][iteration % C_T_len[rotor_idx]] = rotor_state.C_Mx
				average_C_My_arrays[rotor_idx][iteration % C_T_len[rotor_idx]] = rotor_state.C_My
				average_chi_arrays[rotor_idx][iteration % C_T_len[rotor_idx]] = vehicle.inflows[rotor_idx].wake_skew()*180.0/math.pi
				
				average_theta_arrays[rotor_idx][iteration % C_T_len[rotor_idx]] = thetas[rotor_idx,1]
				average_theta_1c_arrays[rotor_idx][iteration % C_T_len[rotor_idx]] = theta_1c[rotor_idx]
				average_theta_1s_arrays[rotor_idx][iteration % C_T_len[rotor_idx]] = theta_1s[rotor_idx]

				average_C_Ts[rotor_idx] = np.sum(average_C_T_arrays[rotor_idx])/C_T_len[rotor_idx]
				average_C_Mxs[rotor_idx] = np.sum(average_C_Mx_arrays[rotor_idx])/C_T_len[rotor_idx]
				average_C_Mys[rotor_idx] = np.sum(average_C_My_arrays[rotor_idx])/C_T_len[rotor_idx]
				average_C_Qs[rotor_idx] = np.sum(average_C_Q_arrays[rotor_idx])/C_T_len[rotor_idx]

				average_theta[rotor_idx] = np.sum(average_theta_arrays[rotor_idx])/C_T_len[rotor_idx]
				average_theta_1c[rotor_idx] = np.sum(average_theta_1c_arrays[rotor_idx])/C_T_len[rotor_idx]
				average_theta_1s[rotor_idx] = np.sum(average_theta_1s_arrays[rotor_idx])/C_T_len[rotor_idx]

				average_Ts[rotor_idx] = average_C_Ts[rotor_idx]*atmo.density*math.pi*vehicle.aircraft.rotors[rotor_idx].radius**4.0*abs(omegas[rotor_idx])**2.0
				average_Mxs[rotor_idx] = average_C_Mxs[rotor_idx]*atmo.density*math.pi*vehicle.aircraft.rotors[rotor_idx].radius**5.0*abs(omegas[rotor_idx])**2.0
				average_Mys[rotor_idx] = -average_C_Mys[rotor_idx]*atmo.density*math.pi*vehicle.aircraft.rotors[rotor_idx].radius**5.0*abs(omegas[rotor_idx])**2.0
				average_Qs[rotor_idx] = average_C_Qs[rotor_idx]*atmo.density*math.pi*vehicle.aircraft.rotors[rotor_idx].radius**5.0*abs(omegas[rotor_idx])**3.0

				average_chis[rotor_idx] = np.sum(average_chi_arrays[rotor_idx])/C_T_len[rotor_idx]

			if moment_trim:
				for rotor_idx, rotor_state in enumerate(vehicle.ac_state.rotor_states):
					curr_c_mxs[rotor_idx] = rotor_state.C_Mx*atmo.density*math.pi*vehicle.aircraft.rotors[rotor_idx].radius**5.0*abs(omegas[rotor_idx])**2.0
					curr_c_mys[rotor_idx] = rotor_state.C_My*atmo.density*math.pi*vehicle.aircraft.rotors[rotor_idx].radius**5.0*abs(omegas[rotor_idx])**2.0
				
				def moment_trim_ode(theta, tau, c_mx_bar, c_my_bar, K, curr_c_mx, curr_c_my, N, rotor_idx):
					a = 8.0/(math.pi*vehicle.aircraft.rotors[rotor_idx].solidity)
					b = flight_condition['moment_trim_const'][rotor_idx]*(1.0/(math.pi*vehicle.aircraft.rotors[rotor_idx].solidity))
					return np.asarray([
						1.0/(tau[0])*(K[0]*(-a*(c_mx_bar - curr_c_mx)/N + b*(c_my_bar - curr_c_my)/N) - theta[0]),
						theta[0],
						-1.0/(tau[1])*(K[1]*(a*(c_my_bar - curr_c_my)/N + b*(c_mx_bar - curr_c_mx)/N) + theta[2]),
						theta[2]
					])

				for rotor_idx, rotor_state in enumerate(vehicle.ac_state.rotor_states):
					N = atmo.density*math.pi*vehicle.aircraft.rotors[rotor_idx].radius**5.0*abs(omegas[rotor_idx])**2.0

					moment_thetas[rotor_idx,:] = moment_thetas[rotor_idx,:] + dt*abs(omegas[rotor_idx])*moment_trim_ode(moment_thetas[rotor_idx,:], moment_taus[rotor_idx,:], c_mx_bars[rotor_idx], c_my_bars[rotor_idx], moment_Ks[rotor_idx,:], curr_c_mxs[rotor_idx], -curr_c_mys[rotor_idx], N, rotor_idx)
					theta_1s[rotor_idx] = moment_thetas[rotor_idx,1]
					theta_1c[rotor_idx] = moment_thetas[rotor_idx,3]

			if trim:
				for rotor_idx, rotor_state in enumerate(vehicle.ac_state.rotor_states):
					curr_c_ts[rotor_idx] = rotor_state.C_T

				def trim_ode(theta, tau, c_t_bar, K, curr_c_t):
					return np.asarray([
						1/tau*(K*6/(2.0*math.pi*vehicle.aircraft.rotors[0].solidity)*(c_t_bar - curr_c_t) - theta[0]),
						theta[0]
					])

				for rotor_idx, rotor_state in enumerate(vehicle.ac_state.rotor_states):
					if trim_algo == 'he':
						thetas[rotor_idx,:] = thetas[rotor_idx,:] + dt*abs(omegas[rotor_idx])*trim_ode(thetas[rotor_idx,:], taus[rotor_idx], c_t_bars[rotor_idx], Ks[rotor_idx], curr_c_ts[rotor_idx])
					elif trim_algo == 'lympany':
						thetas[rotor_idx,0] = thetas[rotor_idx,0] + 3.0/(math.pi*vehicle.aircraft.rotors[rotor_idx].solidity)*(c_t_bars[rotor_idx] - curr_c_ts[rotor_idx])

			for rotor_idx, rotor_state in enumerate(vehicle.ac_state.rotor_states):
				for blade_idx in range(num_blades[rotor_idx]):
					blade_azimuth = vehicle.input_state.rotor_inputs[rotor_idx].azimuth + vehicle.aircraft.rotors[rotor_idx].blades[blade_idx].azimuth_offset
					cos_azimuth = math.cos(blade_azimuth)
					sin_azimuth = math.sin(blade_azimuth)

					sin3_azimuth = math.cos(3.0*(blade_azimuth) - (psi_3 - math.pi))

					if rotor_idx == 0:
						vehicle.input_state.rotor_inputs[rotor_idx].blade_pitches[blade_idx] = thetas[rotor_idx, 1] + theta_1c[rotor_idx]*cos_azimuth + theta_1s[rotor_idx]*sin_azimuth + theta_3*sin3_azimuth
						if elastic_twist:
					  		vehicle.input_state.rotor_inputs[rotor_idx].blade_pitches[blade_idx] = vehicle.input_state.rotor_inputs[rotor_idx].blade_pitches[blade_idx] + elastic_twist(blade_azimuth - math.pi)
					else:
						#vehicle.input_state.rotor_inputs[rotor_idx].blade_pitches[blade_idx] = math.copysign(1, omegas[rotor_idx])*(thetas[rotor_idx, 1] + theta_1c[rotor_idx]*cos_azimuth + theta_1s[rotor_idx]*sin_azimuth + theta_3*sin3_azimuth)
						vehicle.input_state.rotor_inputs[rotor_idx].blade_pitches[blade_idx] = thetas[rotor_idx, 1] + theta_1c[rotor_idx]*cos_azimuth + theta_1s[rotor_idx]*sin_azimuth + theta_3*sin3_azimuth

			if blade_flapping is not None:
				for rotor_idx, rotor_state in enumerate(vehicle.ac_state.rotor_states):
					for blade_idx in range(num_blades[rotor_idx]):
						blade_azimuth = vehicle.input_state.rotor_inputs[rotor_idx].azimuth + vehicle.aircraft.rotors[rotor_idx].blades[blade_idx].azimuth_offset
						(h, h_star) = blade_flapping(blade_azimuth - math.pi)

						vehicle.input_state.rotor_inputs[rotor_idx].blade_flapping_rate[blade_idx] = h_star

			loading_data.time = dt*acoustic_iteration

			if converged:
				if write_wake and (converged_revolutions >= (post_conv_revolutions - 1)):
					for rotor_idx in range(num_rotors):
						write_rotor_vtu(f"{vtu_output_path}/rotor", acoustic_iteration, rotor_idx, vtk_rotors[rotor_idx], vehicle.ac_state.rotor_states[rotor_idx], vehicle.input_state.rotor_inputs[rotor_idx], vehicle.aircraft.rotors[rotor_idx])

						write_wake_vtu(f"{vtu_output_path}/wake", acoustic_iteration, vtk_wake, vehicle.wake_history.history[0])

				if iteration%(post_conv_revolutions*iter_per_rev) == 0 and start_recording:
					done_recording = True

				if iteration%(post_conv_revolutions*iter_per_rev) == 0 and not done_recording:
					start_recording = True

				for rotor_idx, rotor in enumerate(vehicle.ac_state.rotor_states):

					blade_idx = 2   
					# What's the point of making blade_idx 2??
					if rotor_idx > 0:
						blade_idx = 1

					if start_recording and not done_recording:
						for t_idx in range(len(target_span_elements)):
							rotor_idx  # Nitya: WHY?
							dC_L = vehicle.ac_state.rotor_states[rotor_idx].blade_states[blade_idx].chunks[target_span_chunk_index[rotor_idx][t_idx]].dC_L[target_span_element_index[rotor_idx][t_idx]]
							aoa_eff = vehicle.ac_state.rotor_states[rotor_idx].blade_states[blade_idx].chunks[target_span_chunk_index[rotor_idx][t_idx]].aoa_eff[target_span_element_index[rotor_idx][t_idx]]
							aoa = vehicle.ac_state.rotor_states[rotor_idx].blade_states[blade_idx].chunks[target_span_chunk_index[rotor_idx][t_idx]].aoa[target_span_element_index[rotor_idx][t_idx]]
							u_p = vehicle.ac_state.rotor_states[rotor_idx].blade_states[blade_idx].chunks[target_span_chunk_index[rotor_idx][t_idx]].u_p[target_span_element_index[rotor_idx][t_idx]]
							inflow_angle = vehicle.ac_state.rotor_states[rotor_idx].blade_states[blade_idx].chunks[target_span_chunk_index[rotor_idx][t_idx]].inflow_angle[target_span_element_index[rotor_idx][t_idx]]
							theta = vehicle.ac_state.rotor_states[rotor_idx].blade_states[blade_idx].chunks[target_span_chunk_index[rotor_idx][t_idx]].theta[target_span_element_index[rotor_idx][t_idx]]
							gamma = vehicle.ac_state.rotor_states[rotor_idx].blade_states[blade_idx].chunks[target_span_chunk_index[rotor_idx][t_idx]].gamma[target_span_element_index[rotor_idx][t_idx]]

							span_element_loading[rotor_idx, t_idx, spanwise_element_iteration[rotor_idx]] = dC_L*atmo.density*math.pi*vehicle.aircraft.rotors[rotor_idx].radius**3.0*abs(omegas[rotor_idx])**2.0/(0.5*atmo.density*flight_condition["sos"]**2.0*vehicle.aircraft.rotors[rotor_idx].blades[blade_idx].average_chord)
							span_element_aoa_eff[rotor_idx, t_idx, spanwise_element_iteration[rotor_idx]] = aoa_eff*(180.0/math.pi)
							span_element_aoa[rotor_idx, t_idx, spanwise_element_iteration[rotor_idx]] = aoa*(180.0/math.pi)
							span_element_up[rotor_idx, t_idx, spanwise_element_iteration[rotor_idx]] = u_p
							span_element_inflow_angle[rotor_idx, t_idx, spanwise_element_iteration[rotor_idx]] = inflow_angle*(180.0/math.pi)
							span_element_theta[rotor_idx, t_idx, spanwise_element_iteration[rotor_idx]] = theta*(180.0/math.pi)
							span_element_gamma[rotor_idx, t_idx, spanwise_element_iteration[rotor_idx]] = gamma

						spanwise_element_iteration[rotor_idx] = spanwise_element_iteration[rotor_idx] + 1

					for blade_idx, blade in enumerate(rotor.blade_states):

						blade_azimuth = vehicle.input_state.rotor_inputs[rotor_idx].azimuth + vehicle.aircraft.rotors[rotor_idx].blades[blade_idx].azimuth_offset
						
						cos_azimuth = math.cos(blade_azimuth)
						sin_azimuth = math.sin(blade_azimuth)

						sin3_azimuth = math.cos(3.0*(blade_azimuth) - (psi_3 - math.pi))

						collective_pitch_array[rotor_idx, acoustic_iteration] = thetas[rotor_idx, 1]
						sin_pitch_array[rotor_idx, blade_idx, acoustic_iteration] = theta_1s[rotor_idx]*sin_azimuth
						cos_pitch_array[rotor_idx, blade_idx, acoustic_iteration] = theta_1c[rotor_idx]*cos_azimuth
						hhc_pitch_array[rotor_idx, blade_idx, acoustic_iteration] = theta_3*sin3_azimuth

						blade_twist_azimuth[rotor_idx, blade_idx, acoustic_iteration] = acoustic_iteration*d_psi
						if elastic_twist:
							elastic_twist_array[rotor_idx, blade_idx, acoustic_iteration] = elastic_twist(blade_azimuth - math.pi)

						blade_twist_array[rotor_idx, blade_idx, acoustic_iteration] = vehicle.input_state.rotor_inputs[rotor_idx].blade_pitches[blade_idx]

						if blade_flapping is not None:
							(h, h_star) = blade_flapping(blade_azimuth - math.pi)
							blade_flapping_array[rotor_idx, blade_idx, acoustic_iteration] = h
							blade_flapping_der_array[rotor_idx, blade_idx, acoustic_iteration] = h_star

						if track_wake_element:
							for t_idx, target_y_slice in enumerate(target_y_slices):
								fill_wake_y_component(vehicle.wake_history.history[0].rotor_wakes[rotor_idx].tip_vortices[blade_idx], temp_wake_array_y)

								if (abs(temp_wake_array_y[0] - (vehicle.aircraft.rotors[rotor_idx].radius*target_y_slice)) < 0.05) and (not wake_element_found[rotor_idx][t_idx]) and (((blade.azimuth <= 1.02*0.5*math.pi) and (blade.azimuth >= -1.02*0.5*math.pi)) or ((blade.azimuth <= 1.02*5.0/2.0*math.pi) and (blade.azimuth >= 1.02*3.0/2.0*math.pi))):
									num_chunks = blade.chunks.len()
									log_file.write(f"Found wake element to track for slice {target_y_slice} at position: {blade.chunks[num_chunks - 1].x[7]}, {temp_wake_array_y[0]}, {blade.chunks[num_chunks - 1].z[7]}\n")
									wake_element_found[rotor_idx][t_idx] = True
									wake_element_blade[rotor_idx, t_idx] = blade_idx

						fill_dC_Tf(blade, z_loading)
						# fill_dC_cf(blade, y_loading)

						z_loading = -z_loading*atmo.density*math.pi*vehicle.aircraft.rotors[rotor_idx].radius**3.0*abs(omegas[rotor_idx])**2.0
						z_loading = z_loading.astype(dtype=np.single)
						# Nitya: converting z_loading back to single. For some reason, on the lab computer it's getting converted to double!!

						loading_data.set_z_loading_array(z_loading)
						loading_data.set_y_loading_array(y_loading)
						loading_data.set_x_loading_array(x_loading)
						append_loading_data(loading_files[rotor_idx][blade_idx], loading_data)

						# Nitya, 09.16
						for blade_idx2 in range(0,rotor.blade_states.length()):
							wake_idx[blade_idx][blade_idx2][acoustic_iteration].append(get_interactionPt_wake_idx(vehicle.wake_history.history[0].rotor_wakes[0].blade_vortex_interaction[blade_idx].tip_vortex_interaction[blade_idx2]))
							blade_directionVec[blade_idx][blade_idx2][acoustic_iteration].append(get_interaction_point_r_blade(vehicle.wake_history.history[0].rotor_wakes[0].blade_vortex_interaction[blade_idx].tip_vortex_interaction[blade_idx2]))
							vortex_directionVec[blade_idx][blade_idx2][acoustic_iteration].append(get_interaction_point_r_vortex(vehicle.wake_history.history[0].rotor_wakes[0].blade_vortex_interaction[blade_idx].tip_vortex_interaction[blade_idx2]))
							
				if track_wake_element:
					for t_idx in range(len(target_y_slices)):
						for rotor_idx in range(num_rotors):
							if wake_element_found[rotor_idx][t_idx]:
								fill_wake_x_component(vehicle.wake_history.history[0].rotor_wakes[0].tip_vortices[wake_element_blade[rotor_idx, t_idx]], temp_wake_array_x)
								fill_wake_z_component(vehicle.wake_history.history[0].rotor_wakes[0].tip_vortices[wake_element_blade[rotor_idx, t_idx]], temp_wake_array_z)
								fill_wake_r_c_component(vehicle.wake_history.history[0].rotor_wakes[0].tip_vortices[wake_element_blade[rotor_idx, t_idx]], temp_wake_array_r_c)

								if wake_element_index[rotor_idx, t_idx] < min(int(post_conv_revolutions*iter_per_rev), wake_lengths[0]):
									wake_element_trajectory[rotor_idx, t_idx, 0, wake_element_index[rotor_idx, t_idx]] = temp_wake_array_x[wake_element_index[rotor_idx, t_idx]]/vehicle.aircraft.rotors[rotor_idx].radius
									wake_element_trajectory[rotor_idx, t_idx, 1, wake_element_index[rotor_idx, t_idx]] = temp_wake_array_z[wake_element_index[rotor_idx, t_idx]]/vehicle.aircraft.rotors[rotor_idx].radius
									wake_element_core_size[rotor_idx, t_idx, wake_element_index[rotor_idx, t_idx]] = temp_wake_array_r_c[wake_element_index[rotor_idx, t_idx]]

									wake_element_index[rotor_idx, t_idx] = wake_element_index[rotor_idx, t_idx] + 1

				acoustic_iteration = acoustic_iteration + 1

			iteration = iteration + 1
		for rotor_idx in range(num_rotors):
			for blade_idx in range(num_blades[rotor_idx]):
				close_loading_file(loading_files[rotor_idx][blade_idx])

	result_dictionary = {}

	result_dictionary['blade_twist_array'] = blade_twist_array
	result_dictionary['blade_twist_azimuth'] = blade_twist_azimuth

	result_dictionary['collective_pitch_array'] = collective_pitch_array
	result_dictionary['sin_pitch_array'] = sin_pitch_array
	result_dictionary['cos_pitch_array'] = cos_pitch_array
	result_dictionary['hhc_pitch_array'] = hhc_pitch_array

	# Nitya, 09.16
	result_dictionary['wake_idx'] = wake_idx
	#result_dictionary['wake_miss_dist'] = wake_miss_dist
	#result_dictionary['bladeSec_idx'] = bladeSec_idx
	result_dictionary['blade_directionVec'] = blade_directionVec
	result_dictionary['vortex_directionVec'] = vortex_directionVec

	if elastic_twist is not None:
		result_dictionary['elastic_twist_array'] = elastic_twist_array

	if blade_flapping is not None:
		result_dictionary['blade_flapping_array'] = blade_flapping_array
		result_dictionary['blade_flapping_der_array'] = blade_flapping_der_array

	if track_wake_element:
		result_dictionary['wake_element_index'] = wake_element_index
		result_dictionary['target_y_slices'] = target_y_slices
		result_dictionary["wake_element_trajectory"] = wake_element_trajectory
		result_dictionary["wake_element_core_size"] = wake_element_core_size

	if track_span_element:
		result_dictionary['span_element_loading'] = span_element_loading
		result_dictionary['span_element_aoa_eff'] = span_element_aoa_eff
		result_dictionary['span_element_aoa'] = span_element_aoa
		result_dictionary['span_element_up'] = span_element_up
		result_dictionary['span_element_inflow_angle'] = span_element_inflow_angle
		result_dictionary['span_element_theta'] = span_element_theta
		result_dictionary['span_element_gamma'] = span_element_gamma

	namelists = []

	if (acoustics is not None) and (observer is not None):

		min_omega = np.min(np.abs(omegas))

		tau_min = (1.0*math.pi/min_omega)
		tau_max = tau_min + (post_conv_revolutions - 1)*(2.0*math.pi/min_omega)

		min_obs_dist = math.inf
		max_obs_dist = -math.inf

		dist_multiplier = vehicle.aircraft.rotors[observer["reference_rotor"]].radius if observer["radii_relative"] else 1.0

		if observer["type"] == "plane":
			obs_x = np.linspace(observer['min'][0], observer['max'][0], observer['nb'][0])
			obs_y = np.linspace(observer['min'][1], observer['max'][1], observer['nb'][1])
			obs_z = np.linspace(observer['min'][2], observer['max'][2], observer['nb'][2])

			for x, y, z in zip(obs_x, obs_y, obs_z):
				dist = dist_multiplier*math.sqrt(x*x + y*y + z*z)
				min_obs_dist = min(min_obs_dist, dist)
				max_obs_dist = max(max_obs_dist, dist)

		elif observer["type"] == "sphere":
			min_obs_dist = observer["radius"]*dist_multiplier
			max_obs_dist = observer["radius"]*dist_multiplier
		elif observer["type"] == "points":
			for x, y, z in zip(observer['x'], observer['y'], observer['z']):
				x = x*dist_multiplier
				y = y*dist_multiplier
				z = z*dist_multiplier

				for rotor_idx, rotor in enumerate(vehicle.aircraft.rotors):
					for corner_coords in [(rotor.radius, rotor.radius), (rotor.radius, -rotor.radius), (-rotor.radius, rotor.radius), (-rotor.radius, -rotor.radius)]:
						corner = FVec3([
							x - (rotor.frame.global_position()[0] + corner_coords[0]),
							y - (rotor.frame.global_position()[1] + corner_coords[1]),
							z - rotor.frame.global_position()[2]
						])

						dist = math.sqrt(corner[0]*corner[0] + corner[1]*corner[1] + corner[2]*corner[2])

						min_obs_dist = min(min_obs_dist, dist)
						max_obs_dist = max(max_obs_dist, dist)
		elif observer["type"] == "single_point":
			x = observer["x"]*dist_multiplier
			y = observer["y"]*dist_multiplier
			z = observer["z"]*dist_multiplier
			
			for rotor_idx, rotor in enumerate(vehicle.aircraft.rotors):
				for corner_coords in [(rotor.radius, rotor.radius), (rotor.radius, -rotor.radius), (-rotor.radius, rotor.radius), (-rotor.radius, -rotor.radius)]:
					corner = FVec3([
						x - (rotor.frame.global_position()[0] + corner_coords[0]),
						y - (rotor.frame.global_position()[1] + corner_coords[1]),
						z - rotor.frame.global_position()[2]
					])

					dist = math.sqrt(corner[0]*corner[0] + corner[1]*corner[1] + corner[2]*corner[2])

					min_obs_dist = min(min_obs_dist, dist)
					max_obs_dist = max(max_obs_dist, dist)
			
		t_min = tau_min + min_obs_dist/flight_condition["sos"]
		t_max = tau_max + max_obs_dist/flight_condition["sos"]

		nt = int(round((t_max - t_min)/dt))

		for rotor_idx in range(num_rotors):
			wopwop_case_path = f'{wopwop_output_path}/rotor_{rotor_idx}/'

			if not path.isdir(wopwop_case_path):
				makedirs(wopwop_case_path, exist_ok=True)

			namelist = wopwop_input_files_generator.generate_wopwop_namelist(
				atmo,
				dt,
				flight_condition["V_inf"],
				acoustic_iteration,
				flight_condition["aoa"]*(math.pi/180.0),
				t_min,
				t_max,
				nt,
				observer,
				acoustics,
				wopwop_data_path,
				flight_condition["sos"],
				vehicle.aircraft,
				[vehicle.aircraft.rotors[rotor_idx]],
				wopwop_motion,
				vehicle.input_state,
				wopwop_case_path
			)

			namelists.append(namelist)

		wopwop_case_path = f'{wopwop_output_path}/full_system/'

		if not path.isdir(wopwop_case_path):
			makedirs(wopwop_case_path, exist_ok=True)

		namelist = wopwop_input_files_generator.generate_wopwop_namelist(
			atmo,
			dt,
			flight_condition["V_inf"],
			acoustic_iteration,
			flight_condition["aoa"]*(math.pi/180.0),
			t_min,
			t_max,
			nt,
			observer,
			acoustics,
			wopwop_data_path,
			flight_condition["sos"],
			vehicle.aircraft,
			vehicle.aircraft.rotors,
			wopwop_motion,
			vehicle.input_state,
			wopwop_case_path
		)

		namelists.append(namelist)

	return average_C_Ts, namelists, result_dictionary
