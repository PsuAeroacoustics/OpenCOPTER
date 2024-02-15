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

from os import path, makedirs

def flapping_at_azimuth(a, b, w, azimuth):
	h = 0
	h_star = 0

	for idx in range(len(a)):
		cos = math.cos(w*float(idx)*azimuth)
		sin = math.sin(w*float(idx)*azimuth)

		h = h + a[idx]*cos + b[idx - 1]*sin

		if idx > 0:
			h_star = h_star + (-w*float(idx)*a[idx]*sin) + w*float(idx)*b[idx - 1]*cos

	return (h, h_star)

def simulate_aircraft(log_file, vehicle, atmo, r, elements, write_wake, vtu_output_path, wopwop_output_path, do_compute, flight_condition, computational_parameters, observer, acoustics, wake_lengths):

	if not path.isdir(wopwop_output_path):
		makedirs(wopwop_output_path, exist_ok=True)

	if not path.isdir(vtu_output_path):
		makedirs(vtu_output_path, exist_ok=True)

	num_rotors = vehicle.input_state.rotor_inputs.length()

	omegas = np.asarray([vehicle.input_state.rotor_inputs[r_idx].angular_velocity for r_idx in range(num_rotors)])

	log_file.write(f'num_rotors: {num_rotors}\n')
	num_blades = [vehicle.aircraft.rotors[r_idx].blades.length() for r_idx in range(vehicle.aircraft.rotors.length())]

	d_psi = computational_parameters['d_psi']

	dt = d_psi*(math.pi/180.0)/max(abs(omegas))
	iter_per_rev = 360/d_psi

	vtk_rotors = [build_base_vtu_rotor(vehicle.aircraft.rotors[r_idx]) for r_idx in range(num_rotors)]
	vtk_wake = build_base_vtu_wake(vehicle.wake_history.history[0])

	C_T_len = int(round(2.0*math.pi/(dt*max(abs(omegas)))))

	num_rotors = vehicle.input_state.rotor_inputs.length()
	average_C_T_arrays = [np.zeros(C_T_len) for _ in range(num_rotors)]
	average_C_Mx_arrays = [np.zeros(C_T_len) for _ in range(num_rotors)]
	average_C_My_arrays = [np.zeros(C_T_len) for _ in range(num_rotors)]
	average_C_Ts = np.zeros(num_rotors)
	average_C_Mxs = np.zeros(num_rotors)
	average_C_Mys = np.zeros(num_rotors)
	average_Ts = np.zeros(num_rotors)
	average_Mxs = np.zeros(num_rotors)
	average_Mys = np.zeros(num_rotors)

	average_chi_arrays = [np.zeros(C_T_len) for _ in range(num_rotors)]
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
		loading_files = [[wopwop_input_files_generator.build_wopwop_loading(r_idx, blade_idx, int(round(post_conv_revolutions*iter_per_rev)), r, naca0012_xsection, wopwop_data_path) for blade_idx in range(num_blades[r_idx])] for r_idx in range(num_rotors)]

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
		c_t_bars =  [flight_condition["T"][r_idx]/(atmo.density*math.pi*vehicle.aircraft.rotors[r_idx].radius**4.0*abs(omegas[r_idx])**2.0) for r_idx in range(num_rotors)]

	if "c_mx" in flight_condition and "c_my" in flight_condition:
		moment_trim = True
		c_mx_bars = flight_condition["c_mx"]
		c_my_bars = -flight_condition["c_my"]

	elif "Mx" in flight_condition and "My" in flight_condition:
		moment_trim = True
		c_mx_bars =  [flight_condition["Mx"][r_idx] for r_idx in range(num_rotors)]
		c_my_bars =  [-flight_condition["My"][r_idx] for r_idx in range(num_rotors)]

		log_file.write(f'c_mx_bars: {c_mx_bars}, c_my_bars: {c_my_bars}\n')

	if "trim_algo" in computational_parameters:
		trim_algo = computational_parameters["trim_algo"]

	Ks = 0.004*np.ones(num_rotors)
	taus = 12*np.ones(num_rotors)
	moment_Ks = np.ones((num_rotors, 2))
	moment_taus = np.ones((num_rotors, 2))

	# HART II aerodas constants
	moment_Ks[0, 0] = 0.03
	moment_Ks[0, 1] = 0.03

	moment_taus[0, 0] = 0.8
	moment_taus[0, 1] = 0.8

	thetas = np.zeros((num_rotors, 2))
	last_thetas = np.zeros(num_rotors)
	moment_thetas = np.zeros((num_rotors, 2*2))

	theta_1s = np.zeros(num_rotors)
	theta_1c = np.zeros(num_rotors)

	for r_idx in range(num_rotors):
		last_thetas[r_idx] = 10000
		thetas[r_idx, 1] = vehicle.input_state.rotor_inputs[r_idx].blade_pitches[0]

	curr_c_ts = np.zeros(num_rotors)
	curr_c_mxs = np.zeros(num_rotors)
	curr_c_mys = np.zeros(num_rotors)

	last_wake_z = [None for _ in range(num_rotors)]
	wake_l2 = [None for _ in range(num_rotors)]
	sim_done = False
	converged = False
	converged_revolutions = 0
	iteration = 0
	acoustic_iteration = 0
	aoa_update_iter = 0

	start_time = time.perf_counter_ns()

	theta_3 = 0
	psi_3 = 0

	convergence_type = 'wake'
	if 'convergence_type' in computational_parameters:
		convergence_type = computational_parameters['convergence_type']

	if 'psi_3' in flight_condition:
		psi_3 = flight_condition['psi_3']*(math.pi/180.0)
	
	if 'theta_3' in flight_condition:
		theta_3 = flight_condition['theta_3']*(math.pi/180.0)

	if do_compute:

		blade_flapping = None

		if "flapping_coefficients" in flight_condition:
			blade_flapping = lambda az: flapping_at_azimuth(flight_condition["flapping_coefficients"]["cos"], flight_condition["flapping_coefficients"]["sin"], flight_condition["flapping_coefficients"]["w"], az)

		z_loading = np.zeros(elements, dtype=np.single)
		x_loading = np.zeros(elements, dtype=np.single)
		y_loading = np.zeros(elements, dtype=np.single)

		rotor_aoas = [vehicle.input_state.rotor_inputs[r_idx].angle_of_attack for r_idx in range(num_rotors)]
		rotor_aoa = rotor_aoas[0]

		while not sim_done:

			#if iteration % 1 == 0:
			if iteration % iter_per_rev == 0:
				max_l2 = 1000
				if iteration > 0:
					for r_idx in range(num_rotors):
						if convergence_type == 'wake':
							wake_z = np.asarray(get_wake_z_component(vehicle.wake_history.history[0].rotor_wakes[r_idx].tip_vortices[0]))
							wake_l2[r_idx] = np.sqrt(np.mean(np.power((wake_z - last_wake_z[r_idx]), 2.0)))

							max_l2 = np.abs(wake_l2).max()
						elif convergence_type == 'trim':
							max_l2 = np.abs(thetas[r_idx,1] - last_thetas[r_idx])
							last_thetas[r_idx] = thetas[r_idx,1]

					if (max_l2 is not None) and (max_l2 <= computational_parameters["convergence_criteria"]):
						if not converged:
							log_file.write(f"Simulation reached convergence criteria\n")

						converged = True
						aoa_update_iter = aoa_update_iter + 1

						# if convergence_type == 'trim':
						# 	trim = False
							
					else:
						acoustic_iteration = 0
						converged_revolutions = 0
						for rotor_idx, rotor in enumerate(vehicle.ac_state.rotor_states):
							for blade_idx, blade in enumerate(rotor.blade_states):
								restart_loading_file(loading_files[rotor_idx][blade_idx])

						converged = False

				for r_idx in range(num_rotors):
					last_wake_z[r_idx] = np.asarray(get_wake_z_component(vehicle.wake_history.history[0].rotor_wakes[r_idx].tip_vortices[0]))

				now = time.perf_counter_ns()
				elapsed = now - start_time

				log_file.write(
					f'{elapsed/(1000**3):.5f}: rotor rev: {iteration/iter_per_rev:.3f},'
					+''.join([f' C_T{r_idx}: {average_C_Ts[r_idx]:.4f},' for r_idx in range(num_rotors)])
					+''.join([f' T{r_idx}: {average_Ts[r_idx]:.4f},' for r_idx in range(num_rotors)])
					+''.join([f' Mx_{r_idx}: {average_Mxs[r_idx]:.4f},' for r_idx in range(num_rotors)])
					+''.join([f' My_{r_idx}: {average_Mys[r_idx]:.4f},' for r_idx in range(num_rotors)])
					+''.join([f' θ{r_idx}: {thetas[r_idx,1]*(180.0/math.pi):.4f},' for r_idx in range(num_rotors)])
					+''.join([f' θ1c_{r_idx}: {theta_1c[r_idx]*(180.0/math.pi):.4f},' for r_idx in range(num_rotors)])
					+''.join([f' θ1s_{r_idx}: {theta_1s[r_idx]*(180.0/math.pi):.4f},' for r_idx in range(num_rotors)])
					+''.join([f' χ{r_idx}: {average_chis[r_idx]:.2f},' for r_idx in range(num_rotors)])
					+f' combined_C_T: {average_C_Ts.sum():.5f}, max L_2: {max_l2}\n'
				)
				start_time = now
				log_file.flush()


				if converged and not sim_done:
					if converged_revolutions >= post_conv_revolutions:
						sim_done = True

					converged_revolutions = converged_revolutions + 1

			for r_idx in range(vehicle.input_state.rotor_inputs.length()):
				basic_single_rotor_dynamics(vehicle.input_state.rotor_inputs[r_idx], dt)

			step(vehicle.ac_state, vehicle.aircraft, vehicle.input_state, vehicle.inflows, vehicle.wake_history, atmo, iteration, dt)

			for r_idx, rotor_state in enumerate(vehicle.ac_state.rotor_states):
				average_C_T_arrays[r_idx][iteration % C_T_len] = rotor_state.C_T
				average_C_Mx_arrays[r_idx][iteration % C_T_len] = rotor_state.C_Mx
				average_C_My_arrays[r_idx][iteration % C_T_len] = rotor_state.C_My
				average_chi_arrays[r_idx][iteration % C_T_len] = vehicle.inflows[r_idx].wake_skew()*180.0/math.pi

				average_C_Ts[r_idx] = np.sum(average_C_T_arrays[r_idx])/C_T_len
				average_C_Mxs[r_idx] = np.sum(average_C_Mx_arrays[r_idx])/C_T_len
				average_C_Mys[r_idx] = np.sum(average_C_My_arrays[r_idx])/C_T_len

				average_Ts[r_idx] = average_C_Ts[r_idx]*atmo.density*math.pi*vehicle.aircraft.rotors[r_idx].radius**4.0*abs(omegas[r_idx])**2.0
				average_Mxs[r_idx] = average_C_Mxs[r_idx]*atmo.density*math.pi*vehicle.aircraft.rotors[r_idx].radius**5.0*abs(omegas[r_idx])**2.0
				average_Mys[r_idx] = average_C_Mys[r_idx]*atmo.density*math.pi*vehicle.aircraft.rotors[r_idx].radius**5.0*abs(omegas[r_idx])**2.0

				average_chis[r_idx] = np.sum(average_chi_arrays[r_idx])/C_T_len

			if moment_trim:
				for r_idx, rotor_state in enumerate(vehicle.ac_state.rotor_states):
					curr_c_mxs[r_idx] = rotor_state.C_Mx*atmo.density*math.pi*vehicle.aircraft.rotors[r_idx].radius**5.0*abs(omegas[r_idx])**2.0
					curr_c_mys[r_idx] = -rotor_state.C_My*atmo.density*math.pi*vehicle.aircraft.rotors[r_idx].radius**5.0*abs(omegas[r_idx])**2.0
				
				def moment_trim_ode(theta, tau, c_mx_bar, c_my_bar, K, curr_c_mx, curr_c_my, N):
					a = 8.0/(math.pi*vehicle.aircraft.rotors[0].solidity)
					b = flight_condition['moment_trim_const']*(1.0/(math.pi*vehicle.aircraft.rotors[0].solidity))
					return np.asarray([
						1.0/(tau[0])*(K[0]*(-a*(c_mx_bar - curr_c_mx)/N + b*(c_my_bar - curr_c_my)/N) - theta[0]),
						theta[0],
						-1.0/(tau[1])*(K[1]*(a*(c_my_bar - curr_c_my)/N + b*(c_mx_bar - curr_c_mx)/N) + theta[2]),
						theta[2]
					])

				for r_idx, rotor_state in enumerate(vehicle.ac_state.rotor_states):
					N = atmo.density*math.pi*vehicle.aircraft.rotors[r_idx].radius**5.0*abs(omegas[r_idx])**2.0

					moment_thetas[r_idx,:] = moment_thetas[r_idx,:] + moment_trim_ode(moment_thetas[r_idx,:], moment_taus[r_idx,:], c_mx_bars[r_idx], c_my_bars[r_idx], moment_Ks[r_idx,:], curr_c_mxs[r_idx], curr_c_mys[r_idx], N)
					theta_1s[r_idx] = moment_thetas[r_idx,1]
					theta_1c[r_idx] = moment_thetas[r_idx,3]

			if trim:
				for r_idx, rotor_state in enumerate(vehicle.ac_state.rotor_states):
					curr_c_ts[r_idx] = rotor_state.C_T
					

				def trim_ode(theta, tau, c_t_bar, K, curr_c_t):
					return np.asarray([
						1/tau*(K*6/(2.0*math.pi*vehicle.aircraft.rotors[0].solidity)*(c_t_bar - curr_c_t) - theta[0]),
						theta[0]
					])

				for r_idx, rotor_state in enumerate(vehicle.ac_state.rotor_states):
					if trim_algo == 'he':
						thetas[r_idx,:] = thetas[r_idx,:] + trim_ode(thetas[r_idx,:], taus[r_idx], c_t_bars[r_idx], Ks[r_idx], curr_c_ts[r_idx])
					elif trim_algo == 'lympany':
						thetas[r_idx,0] = thetas[r_idx,0] + 3.0/(math.pi*vehicle.aircraft.rotors[r_idx].solidity)*(c_t_bars[r_idx] - curr_c_ts[r_idx])
					
					for b_idx in range(num_blades[r_idx]):
						blade_azimuth = math.fmod(vehicle.input_state.rotor_inputs[r_idx].azimuth + vehicle.aircraft.rotors[r_idx].blades[b_idx].azimuth_offset, 2.0*math.pi)
						cos_azimuth = math.cos(blade_azimuth)
						sin_azimuth = math.sin(blade_azimuth)

						sin3_azimuth = math.sin(3.0*blade_azimuth + psi_3)

						vehicle.input_state.rotor_inputs[r_idx].blade_pitches[b_idx] = thetas[r_idx, 1] + theta_1c[r_idx]*cos_azimuth + theta_1s[r_idx]*sin_azimuth + theta_3*sin3_azimuth

			if blade_flapping is not None:
				for r_idx, rotor_state in enumerate(vehicle.ac_state.rotor_states):
					for b_idx in range(num_blades[r_idx]):
						blade_azimuth = vehicle.input_state.rotor_inputs[r_idx].azimuth + vehicle.aircraft.rotors[r_idx].blades[b_idx].azimuth_offset
						(h, h_star) = blade_flapping(blade_azimuth)

						vehicle.input_state.rotor_inputs[r_idx].blade_flapping[b_idx] = -h - flight_condition["flapping_coefficients"]["coning"]*math.pi/180.0
						vehicle.input_state.rotor_inputs[r_idx].blade_flapping_rate[b_idx] = h_star
						
			loading_data.time = dt*acoustic_iteration

			if converged and (vehicle.input_state.rotor_inputs[r_idx].angle_of_attack == rotor_aoa):
				if write_wake and (converged_revolutions >= (post_conv_revolutions - 1)):
					for r_idx in range(num_rotors):
						write_rotor_vtu(f"{vtu_output_path}/rotor", acoustic_iteration, r_idx, vtk_rotors[r_idx], vehicle.ac_state.rotor_states[r_idx], vehicle.input_state.rotor_inputs[r_idx])

						write_wake_vtu(f"{vtu_output_path}/wake", acoustic_iteration, vtk_wake, vehicle.wake_history.history[0])

				for rotor_idx, rotor in enumerate(vehicle.ac_state.rotor_states):
					for blade_idx, blade in enumerate(rotor.blade_states):
						fill_dC_Nf(blade, z_loading)
						fill_dC_cf(blade, y_loading)
						#z_loading = [-z_load*atmo.density*math.pi*vehicle.aircraft.rotors[rotor_idx].radius**4.0*abs(omegas[rotor_idx])**2.0/vehicle.aircraft.rotors[rotor_idx].blades[blade_idx].average_chord for z_load in z_loading]

						#z_loading = -1.0*np.abs(z_loading)*atmo.density*math.pi*vehicle.aircraft.rotors[rotor_idx].radius**4.0*abs(omegas[rotor_idx])**2.0
						z_loading = -z_loading*atmo.density*math.pi*vehicle.aircraft.rotors[rotor_idx].radius**3.0*abs(omegas[rotor_idx])**2.0
						#y_loading = np.sign(omegas[rotor_idx])*y_loading*atmo.density*math.pi*vehicle.aircraft.rotors[rotor_idx].radius**2.0*abs(omegas[rotor_idx])**2.0

						loading_data.set_z_loading_array(z_loading)
						loading_data.set_y_loading_array(y_loading)
						loading_data.set_x_loading_array(x_loading)

						append_loading_data(loading_files[rotor_idx][blade_idx], loading_data)

				acoustic_iteration = acoustic_iteration + 1

			iteration = iteration + 1
		for rotor_idx in range(num_rotors):
			for blade_idx in range(num_blades[rotor_idx]):
				close_loading_file(loading_files[rotor_idx][blade_idx])

	namelists = []
	collectives = [vehicle.input_state.rotor_inputs[r_idx].blade_pitches[0] for r_idx in range(num_rotors)]

	if (acoustics is not None) and (observer is not None):

		min_omega = np.min(np.abs(omegas))

		tau_min = (1.0*math.pi/min_omega)
		tau_max = tau_min + (post_conv_revolutions - 1)*(2.0*math.pi/min_omega)

		min_obs_dist = math.inf
		max_obs_dist = -math.inf

		max_R = max([vehicle.aircraft.rotors[r_idx].radius for r_idx in range(num_rotors)])
		
		dist_multiplier = max_R if observer["radii_relative"] else 1.0

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
				for r_idx, rotor in enumerate(vehicle.aircraft.rotors):
					for corner_coords in [(max_R, max_R), (max_R, -max_R), (-max_R, max_R), (-max_R, -max_R)]:
						corner = FVec3([
							x - (rotor.origin[0] + corner_coords[0]),
							y - (rotor.origin[1] + corner_coords[1]),
							z - rotor.origin[2]
						])

						dist = dist_multiplier*math.sqrt(corner[0]*corner[0] + corner[1]*corner[1] + corner[2]*corner[2])

						min_obs_dist = min(min_obs_dist, dist)
						max_obs_dist = max(max_obs_dist, dist)

		t_min = tau_min + min_obs_dist/flight_condition["sos"]
		t_max = tau_max + max_obs_dist/flight_condition["sos"]

		nt = int(round((t_max - t_min)/dt))

		for r_idx in range(num_rotors):
			wopwop_case_path = f'{wopwop_output_path}/rotor_{r_idx}/'

			if not path.isdir(wopwop_case_path):
				makedirs(wopwop_case_path, exist_ok=True)

			namelist = wopwop_input_files_generator.generate_wopwop_namelist(
				[vehicle.aircraft.rotors[r_idx].radius for r_idx in range(num_rotors)],
				[vehicle.aircraft.rotors[r_idx].origin for r_idx in range(num_rotors)],
				atmo,
				1,
				num_blades[r_idx],
				[omegas[r_idx]],
				dt,
				flight_condition["V_inf"],
				3*iter_per_rev,
				vehicle.input_state.rotor_inputs[r_idx].angle_of_attack,
				r_idx,
				t_min,
				t_max,
				nt,
				observer,
				acoustics,
				wopwop_data_path,
				flight_condition["sos"],
				[collectives[r_idx]],
				False
			)

			namelists.append(namelist)

		wopwop_case_path = f'{wopwop_output_path}/full_system/'

		if not path.isdir(wopwop_case_path):
			makedirs(wopwop_case_path, exist_ok=True)

		namelist = wopwop_input_files_generator.generate_wopwop_namelist(
			[vehicle.aircraft.rotors[r_idx].radius for r_idx in range(num_rotors)],
			[vehicle.aircraft.rotors[r_idx].origin for r_idx in range(num_rotors)],
			atmo,
			num_rotors,
			num_blades,
			omegas,
			dt,
			flight_condition["V_inf"],
			3*iter_per_rev,
			vehicle.input_state.rotor_inputs[0].angle_of_attack,
			r_idx,
			t_min,
			t_max,
			nt,
			observer,
			acoustics,
			wopwop_data_path,
			flight_condition["sos"],
			collectives,
			True
		)

		namelists.append(namelist)

	return average_C_Ts, namelists
