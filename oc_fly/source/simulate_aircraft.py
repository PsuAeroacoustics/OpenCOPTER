from ast import operator
import sys
import os

sys.path.append(f'{os.path.dirname(os.path.realpath(__file__))}/dependencies/OpenCOPTER')
sys.path.append(f'{os.path.dirname(os.path.realpath(__file__))}/dependencies/wopwopd')

import scipy.optimize as scp
from libopencopter import *
from libwopwopd import *
import wopwop_input_files_generator
import numpy as np
import math
import time
import argparse
from collections.abc import Iterable
import functools

from os import walk, path, makedirs

def simulate_aircraft(log_file, aircraft, ac_state, input_state, inflows, wake_history, atmo, omegas, V_inf, r, twist, AR, d_psi, elements, wopwop_rotor_idx, write_wake, vtu_output_path, wopwop_output_path, do_compute, flight_condition, convergence_criteria):

	if not path.isdir(wopwop_output_path):
		makedirs(wopwop_output_path, exist_ok=True)

	if not path.isdir(vtu_output_path):
		makedirs(vtu_output_path, exist_ok=True)

	num_rotors = input_state.rotor_inputs.length()

	log_file.write(f"wopwop_rotor_indx = {wopwop_rotor_idx}\n")

	if not(0 <= wopwop_rotor_idx <= num_rotors):
		log_file.write(f"Error: wopwop_rotor_indx is {wopwop_rotor_idx} which is more than the number of rotors ({num_rotors}). Check all the inputs. \t terminating\n")
		sys.exit()

	omegas = np.asarray(omegas)
	log_file.write(f'num_rotors: {num_rotors}\n')
	num_blades = aircraft.rotors[0].blades.length()

	dt = d_psi*(math.pi/180.0)/max(abs(omegas))
	iter_per_rev = 360/d_psi

	vtk_rotors = [build_base_vtu_rotor(aircraft.rotors[r_idx]) for r_idx in range(num_rotors)]
	vtk_wake = build_base_vtu_wake(wake_history.history[0])

	C_T_len = int(round(2.0*math.pi/(dt*max(abs(omegas)))))

	num_rotors = input_state.rotor_inputs.length()
	average_C_T_arrays = [np.zeros(C_T_len) for _ in range(num_rotors)]
	average_C_Ts = np.zeros(num_rotors)

	average_chi_arrays = [np.zeros(C_T_len) for _ in range(num_rotors)]
	average_chis = np.zeros(num_rotors)

	loading_data = ZoneVectorData(elements)

	loading_data.set_x_loading_array(np.zeros(elements, dtype=np.single))
	loading_data.set_y_loading_array(np.zeros(elements, dtype=np.single))

	#iterations = int(round(2*math.pi/d_psi*rotor_revs))
	#iterations = int(round(iter_per_rev*rotor_revs))
	#log_file.write("Number of iterations = ", iterations)

	naca0012_xsection = naca0012()

	wopwop_data_path = f'{wopwop_output_path}/data'

	if not path.isdir(wopwop_data_path):
		makedirs(wopwop_data_path, exist_ok=True)

	wopwop_input_files_generator.write_wopwop_geometry(naca0012_xsection, r, twist, aircraft.rotors[wopwop_rotor_idx].radius, AR, wopwop_data_path)

	loading_files = []

	if do_compute:
		loading_files = [[wopwop_input_files_generator.build_wopwop_loading(r_idx, blade_idx, int(round(2*iter_per_rev)), r, naca0012_xsection, wopwop_data_path) for blade_idx in range(num_blades)] for r_idx in range(num_rotors)]

	log_file.write("Performing acoustic resolution simulation\n")
	start_time = time.perf_counter()

	trim = False
	c_t_bars = None
	if "c_t" in flight_condition:
		trim = True
		c_t_bars = flight_condition["c_t"]

	Ks = 0.3*np.ones(num_rotors)
	taus = 2.88*np.ones(num_rotors)

	thetas = np.zeros((num_rotors, 2))

	for r_idx in range(num_rotors):
		thetas[r_idx, 0] = input_state.rotor_inputs[r_idx].blade_pitches[0]
		thetas[r_idx, 1] = 0

	curr_c_ts = np.zeros(2)

	last_wake_z = [None for _ in range(num_rotors)]
	wake_l2 = [None for _ in range(num_rotors)]
	sim_done = False
	converged = False
	converged_revolutions = 0
	iteration = 0
	acoustic_iteration = 0

	if do_compute:
		# x_inflow = np.linspace(-13, 2, 2048)

		# y_slice = np.zeros(chunk_size())
		# x_e_slice = np.zeros(chunk_size())
		# z_slice = np.zeros(chunk_size())

		# x_inflow = x_inflow.reshape((int(2048/chunk_size()), chunk_size()))
		# inflow_induced_velocities = np.zeros((int(2048/chunk_size()), chunk_size()))
		# wake_induced_velocities = np.zeros((int(2048/chunk_size()), chunk_size()))

		rotor_aoas = [input_state.rotor_inputs[r_idx].angle_of_attack for r_idx in range(num_rotors)]

		for r_idx, rotor_aoa in enumerate(rotor_aoas):
			if rotor_aoa < 0.0:
				input_state.rotor_inputs[r_idx].angle_of_attack = 0

		while not sim_done:

			for r_idx, rotor_aoa in enumerate(rotor_aoas):
				if rotor_aoa < 0.0:
					input_state.rotor_inputs[r_idx].angle_of_attack = np.sign(rotor_aoa)*(0.5/360.0)*(math.pi/180.0)*float(iteration)

					if input_state.rotor_inputs[r_idx].angle_of_attack < rotor_aoa:
						input_state.rotor_inputs[r_idx].angle_of_attack = rotor_aoa

			if iteration % iter_per_rev == 0:
				max_l2 = None
				if iteration > 0:
					for r_idx in range(num_rotors):
						wake_z = np.asarray(get_wake_z_component(wake_history.history[0].rotor_wakes[r_idx].tip_vortices[0]))
						wake_l2[r_idx] = np.sqrt(np.mean(np.power((wake_z - last_wake_z[r_idx]), 2.0)))

					max_l2 = np.abs(wake_l2).max()

					if (max_l2 is not None) and (max_l2 <= convergence_criteria):
						converged = True

				for r_idx in range(num_rotors):
					last_wake_z[r_idx] = np.asarray(get_wake_z_component(wake_history.history[0].rotor_wakes[r_idx].tip_vortices[0]))

				now = time.perf_counter()
				elapsed = now - start_time
				#log_file.write(f'{elapsed:.5f}: rotor rev: {iteration/360}, C_T1: {average_C_T1:.8f}, C_T2: {average_C_T2:.8f}, χ1: {average_chi1:.4}, χ2: {average_chi2:.4}')
				log_file.write(f'{elapsed:.5f}: rotor rev: {iteration/iter_per_rev},'+''.join([f' C_T{r_idx}: {average_C_Ts[r_idx]:.8f},' for r_idx in range(num_rotors)])+''.join([f' χ{r_idx}: {average_chis[r_idx]:.4f},' for r_idx in range(num_rotors)])+f' combined_C_T: {average_C_Ts.sum()}, max L_2: {max_l2}\n')
				start_time = now
				log_file.flush()


				if converged and not sim_done:
					if converged_revolutions >= 2:
						sim_done = True

					converged_revolutions = converged_revolutions + 1


			for r_idx in range(input_state.rotor_inputs.length()):
				basic_single_rotor_dynamics(input_state.rotor_inputs[r_idx], dt)

			step(ac_state, aircraft, input_state, inflows, wake_history, atmo, iteration, dt)

			for r_idx, rotor_state in enumerate(ac_state.rotor_states):
				average_C_T_arrays[r_idx][iteration % C_T_len] = rotor_state.C_T
				average_chi_arrays[r_idx][iteration % C_T_len] = inflows[r_idx].wake_skew()*180.0/math.pi

				average_C_Ts[r_idx] = np.sum(average_C_T_arrays[r_idx])/C_T_len
				average_chis[r_idx] = np.sum(average_chi_arrays[r_idx])/C_T_len


			if trim:
				for r_idx, rotor_state in enumerate(ac_state.rotor_states):
					curr_c_ts[r_idx] = rotor_state.C_T

				def trim_ode(theta, tau, c_t_bar, K, curr_c_t):
					return np.asarray([
						1/tau*(K*6/(2.0*math.pi*aircraft.rotors[0].solidity)*(c_t_bar - curr_c_t) - theta[1]),
						theta[1]
					])

				for r_idx, rotor_state in enumerate(ac_state.rotor_states):
					thetas[r_idx,:] = thetas[r_idx,:] + trim_ode(thetas[r_idx,:], taus[r_idx], c_t_bars[r_idx], Ks[r_idx], curr_c_ts[r_idx])
					
					for b_idx in range(num_blades):
						input_state.rotor_inputs[r_idx].blade_pitches[b_idx] = thetas[r_idx, 0]

			loading_data.time = dt*acoustic_iteration

			if converged:

				if write_wake:
					for r_idx in range(num_rotors):
						write_rotor_vtu(f"{vtu_output_path}/rotor", acoustic_iteration, r_idx, vtk_rotors[r_idx], ac_state.rotor_states[r_idx], input_state.rotor_inputs[r_idx])

						write_wake_vtu(f"{vtu_output_path}/wake", acoustic_iteration, vtk_wake, wake_history.history[0])

				# for c_idx, x_chunk in enumerate(x_inflow):
				# 	induced_velocities = compute_wake_induced_velocities(wake_history.history[0], x_chunk, y_slice, z_slice, ac_state, omegas[0], 0, True)

				# 	inflow_induced_velocities[c_idx, :] = inflow_induced_velocities[c_idx, :] + np.asarray(inflows[0].inflow_at(x_chunk, y_slice, z_slice, x_e_slice, input_state.rotor_inputs[0].angle_of_attack))/(2*iter_per_rev)
				# 	wake_induced_velocities[c_idx, :] = wake_induced_velocities[c_idx, :] + np.asarray(induced_velocities.v_z)/(2*iter_per_rev)
	
				for rotor_idx, rotor in enumerate(ac_state.rotor_states):
					for blade_idx, blade in enumerate(rotor.blade_states):
						z_loading = get_dC_T(blade)
						z_loading = [0.5*z_load*atmo.density*math.pi*aircraft.rotors[rotor_idx].radius**4.0*abs(omegas[rotor_idx])**2.0 for z_load in z_loading]
						loading_data.set_z_loading_array(z_loading)

						append_loading_data(loading_files[rotor_idx][blade_idx], loading_data)

				acoustic_iteration = acoustic_iteration + 1

			iteration = iteration + 1
		for rotor_idx in range(num_rotors):
			for blade_idx in range(num_blades):
				close_loading_file(loading_files[rotor_idx][blade_idx])

	acoustic_source_start_iter = 0 #acoustic_iteration - 2*iter_per_rev
	tau_min = 0#acoustic_source_start_iter*dt

	t_min = tau_min + 20*aircraft.rotors[wopwop_rotor_idx].radius/343
	t_max = t_min + (2.0*math.pi/abs(omegas[wopwop_rotor_idx]))

	nt = int(round((t_max - t_min)/dt))

	namelists = []

	for r_idx in range(num_rotors):
		wopwop_case_path = f'{wopwop_output_path}/rotor_{r_idx}/'

		if not path.isdir(wopwop_case_path):
				makedirs(wopwop_case_path, exist_ok=True)

		namelist = wopwop_input_files_generator.generate_wopwop_namelist(
			[aircraft.rotors[r_idx].radius for r_idx in range(num_rotors)],
			[aircraft.rotors[r_idx].origin for r_idx in range(num_rotors)],
			atmo,
			1,
			num_blades,
			omegas,
			dt,
			V_inf,
			acoustic_iteration,
			input_state.rotor_inputs[r_idx].angle_of_attack,
			wopwop_data_path,
			r_idx,
			t_min,
			t_max,
			nt,
			False
		)

		namelists.append(namelist)

	wopwop_case_path = f'{wopwop_output_path}/full_system/'

	if not path.isdir(wopwop_case_path):
		makedirs(wopwop_case_path, exist_ok=True)

	namelist = wopwop_input_files_generator.generate_wopwop_namelist(
		[aircraft.rotors[r_idx].radius for r_idx in range(num_rotors)],
		[aircraft.rotors[r_idx].origin for r_idx in range(num_rotors)],
		atmo,
		num_rotors,
		num_blades,
		omegas,
		dt,
		V_inf,
		acoustic_iteration,
		input_state.rotor_inputs[0].angle_of_attack,
		wopwop_data_path,
		r_idx,
		t_min,
		t_max,
		nt,
		True
	)

	namelists.append(namelist)

	#inflow_induced_velocities = inflow_induced_velocities.reshape(2048)
	#wake_induced_velocities = wake_induced_velocities.reshape(2048)
	#x_inflow = x_inflow.reshape(2048)

	return average_C_Ts, namelists#, x_inflow, inflow_induced_velocities, wake_induced_velocities
