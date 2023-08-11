from ast import operator
import sys
import os

sys.path.append(f'{os.path.dirname(os.path.realpath(__file__))}/dependencies/OpenCOPTER')
sys.path.append(f'{os.path.dirname(os.path.realpath(__file__))}/dependencies/wopwopd')

from libopencopter import *
from libwopwopd import *
import wopwop_input_files_generator
import numpy as np
import math
import time

from os import path, makedirs

def simulate_aircraft(log_file, vehicle, atmo, r, d_psi, elements, write_wake, vtu_output_path, wopwop_output_path, do_compute, flight_condition, convergence_criteria, observer, acoustics):

	if not path.isdir(wopwop_output_path):
		makedirs(wopwop_output_path, exist_ok=True)

	if not path.isdir(vtu_output_path):
		makedirs(vtu_output_path, exist_ok=True)

	num_rotors = vehicle.input_state.rotor_inputs.length()

	omegas = np.asarray([vehicle.input_state.rotor_inputs[r_idx].angular_velocity for r_idx in range(num_rotors)])

	log_file.write(f'num_rotors: {num_rotors}\n')
	num_blades = [vehicle.aircraft.rotors[r_idx].blades.length() for r_idx in range(vehicle.aircraft.rotors.length())]

	dt = d_psi*(math.pi/180.0)/max(abs(omegas))
	iter_per_rev = 360/d_psi

	vtk_rotors = [build_base_vtu_rotor(vehicle.aircraft.rotors[r_idx]) for r_idx in range(num_rotors)]
	vtk_wake = build_base_vtu_wake(vehicle.wake_history.history[0])

	C_T_len = int(round(2.0*math.pi/(dt*max(abs(omegas)))))

	num_rotors = vehicle.input_state.rotor_inputs.length()
	average_C_T_arrays = [np.zeros(C_T_len) for _ in range(num_rotors)]
	average_C_Ts = np.zeros(num_rotors)

	average_chi_arrays = [np.zeros(C_T_len) for _ in range(num_rotors)]
	average_chis = np.zeros(num_rotors)

	loading_data = ZoneVectorData(elements)

	loading_data.set_x_loading_array(np.zeros(elements, dtype=np.single))
	loading_data.set_y_loading_array(np.zeros(elements, dtype=np.single))

	naca0012_xsection = naca0012()

	wopwop_data_path = f'{wopwop_output_path}/data'

	loading_files = []

	if do_compute:
		loading_files = [[wopwop_input_files_generator.build_wopwop_loading(r_idx, blade_idx, int(round(2*iter_per_rev)), r, naca0012_xsection, wopwop_data_path) for blade_idx in range(num_blades[r_idx])] for r_idx in range(num_rotors)]

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
		thetas[r_idx, 0] = vehicle.input_state.rotor_inputs[r_idx].blade_pitches[0]
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

		rotor_aoas = [vehicle.input_state.rotor_inputs[r_idx].angle_of_attack for r_idx in range(num_rotors)]

		for r_idx, rotor_aoa in enumerate(rotor_aoas):
			if rotor_aoa < 0.0:
				vehicle.input_state.rotor_inputs[r_idx].angle_of_attack = 0

		while not sim_done:

			for r_idx, rotor_aoa in enumerate(rotor_aoas):
				if rotor_aoa < 0.0:
					vehicle.input_state.rotor_inputs[r_idx].angle_of_attack = np.sign(rotor_aoa)*(0.5/360.0)*(math.pi/180.0)*float(iteration)

					if vehicle.input_state.rotor_inputs[r_idx].angle_of_attack < rotor_aoa:
						vehicle.input_state.rotor_inputs[r_idx].angle_of_attack = rotor_aoa

			if iteration % iter_per_rev == 0:
				max_l2 = None
				if iteration > 0:
					for r_idx in range(num_rotors):
						wake_z = np.asarray(get_wake_z_component(vehicle.wake_history.history[0].rotor_wakes[r_idx].tip_vortices[0]))
						wake_l2[r_idx] = np.sqrt(np.mean(np.power((wake_z - last_wake_z[r_idx]), 2.0)))

					max_l2 = np.abs(wake_l2).max()

					if (max_l2 is not None) and (max_l2 <= convergence_criteria):
						converged = True

				for r_idx in range(num_rotors):
					last_wake_z[r_idx] = np.asarray(get_wake_z_component(vehicle.wake_history.history[0].rotor_wakes[r_idx].tip_vortices[0]))

				now = time.perf_counter()
				elapsed = now - start_time

				log_file.write(f'{elapsed:.5f}: rotor rev: {iteration/iter_per_rev},'+''.join([f' C_T{r_idx}: {average_C_Ts[r_idx]:.8f},' for r_idx in range(num_rotors)])+''.join([f' χ{r_idx}: {average_chis[r_idx]:.4f},' for r_idx in range(num_rotors)])+f' combined_C_T: {average_C_Ts.sum()}, max L_2: {max_l2}\n')
				start_time = now
				log_file.flush()


				if converged and not sim_done:
					if converged_revolutions >= 2:
						sim_done = True

					converged_revolutions = converged_revolutions + 1


			for r_idx in range(vehicle.input_state.rotor_inputs.length()):
				basic_single_rotor_dynamics(vehicle.input_state.rotor_inputs[r_idx], dt)

			step(vehicle.ac_state, vehicle.aircraft, vehicle.input_state, vehicle.inflows, vehicle.wake_history, atmo, iteration, dt)

			for r_idx, rotor_state in enumerate(vehicle.ac_state.rotor_states):
				average_C_T_arrays[r_idx][iteration % C_T_len] = rotor_state.C_T
				average_chi_arrays[r_idx][iteration % C_T_len] = vehicle.inflows[r_idx].wake_skew()*180.0/math.pi

				average_C_Ts[r_idx] = np.sum(average_C_T_arrays[r_idx])/C_T_len
				average_chis[r_idx] = np.sum(average_chi_arrays[r_idx])/C_T_len

			if trim:
				for r_idx, rotor_state in enumerate(vehicle.ac_state.rotor_states):
					curr_c_ts[r_idx] = rotor_state.C_T

				def trim_ode(theta, tau, c_t_bar, K, curr_c_t):
					return np.asarray([
						1/tau*(K*6/(2.0*math.pi*vehicle.aircraft.rotors[0].solidity)*(c_t_bar - curr_c_t) - theta[1]),
						theta[1]
					])

				for r_idx, rotor_state in enumerate(vehicle.ac_state.rotor_states):
					thetas[r_idx,:] = thetas[r_idx,:] + trim_ode(thetas[r_idx,:], taus[r_idx], c_t_bars[r_idx], Ks[r_idx], curr_c_ts[r_idx])
					
					for b_idx in range(num_blades[r_idx]):
						vehicle.input_state.rotor_inputs[r_idx].blade_pitches[b_idx] = thetas[r_idx, 0]

			loading_data.time = dt*acoustic_iteration

			# if write_wake:
			# 	for r_idx in range(num_rotors):
			# 		write_rotor_vtu(f"{vtu_output_path}/rotor", iteration, r_idx, vtk_rotors[r_idx], ac_state.rotor_states[r_idx], input_state.rotor_inputs[r_idx])

			# 		write_wake_vtu(f"{vtu_output_path}/wake", iteration, vtk_wake, wake_history.history[0])

			if converged:
				if write_wake:
					for r_idx in range(num_rotors):
						write_rotor_vtu(f"{vtu_output_path}/rotor", acoustic_iteration, r_idx, vtk_rotors[r_idx], vehicle.ac_state.rotor_states[r_idx], vehicle.input_state.rotor_inputs[r_idx])

						write_wake_vtu(f"{vtu_output_path}/wake", acoustic_iteration, vtk_wake, vehicle.wake_history.history[0])

				for rotor_idx, rotor in enumerate(vehicle.ac_state.rotor_states):
					for blade_idx, blade in enumerate(rotor.blade_states):
						z_loading = get_dC_T(blade)
						z_loading = [0.5*z_load*atmo.density*math.pi*vehicle.aircraft.rotors[rotor_idx].radius**4.0*abs(omegas[rotor_idx])**2.0 for z_load in z_loading]
						loading_data.set_z_loading_array(z_loading)

						append_loading_data(loading_files[rotor_idx][blade_idx], loading_data)

				acoustic_iteration = acoustic_iteration + 1

			iteration = iteration + 1
		for rotor_idx in range(num_rotors):
			for blade_idx in range(num_blades[rotor_idx]):
				close_loading_file(loading_files[rotor_idx][blade_idx])

	tau_min = 0

	max_R = max([vehicle.aircraft.rotors[r_idx].radius for r_idx in range(num_rotors)])
	min_omega = np.min(np.abs(omegas))

	t_min = tau_min + 20*max_R/343
	t_max = t_min + (2.0*math.pi/min_omega)

	nt = int(round((t_max - t_min)/dt))

	namelists = []

	if (acoustics is not None) and (observer is not None):
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
				omegas,
				dt,
				flight_condition["V_inf"],
				acoustic_iteration,
				vehicle.input_state.rotor_inputs[r_idx].angle_of_attack,
				r_idx,
				t_min,
				t_max,
				nt,
				observer,
				acoustics,
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
			acoustic_iteration,
			vehicle.input_state.rotor_inputs[0].angle_of_attack,
			r_idx,
			t_min,
			t_max,
			nt,
			observer,
			acoustics,
			True
		)

		namelists.append(namelist)

	return average_C_Ts, namelists
