from ast import operator
import sys
import os
#from memory_profiler import profile

sys.path.append(f'{os.path.dirname(os.path.realpath(__file__))}/dependencies/OpenCOPTER')
sys.path.append(f'{os.path.dirname(os.path.realpath(__file__))}/dependencies/wopwopd')

from libopencopter import *
from libwopwopd import *
import wopwop_input_files_generator
import numpy as np
from scipy.integrate import simpson
import scipy.io
import math
import time

from simulated_vehicle import SimulatedVehicle

from os import path, makedirs

TRIM_MODE_COLLECTIVE = 0
TRIM_MODE_SHARED_COLLECTIVE = 1
TRIM_MODE_RPM = 2

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

	for idx in range(len(a)):
		cos = math.cos(float(idx)*azimuth)
		sin = math.sin(float(idx)*azimuth)

		h = h + a[idx]*cos + b[idx - 1]*sin

	return h

def simulate_aircraft(log_file, vehicle: SimulatedVehicle, atmo, elements, write_wake, output_base, vtu_output_path, wopwop_output_path, do_compute, flight_condition, computational_parameters, observer, acoustics, wake_lengths, results, wopwop_motion, geom_directory):
	if not path.isdir(wopwop_output_path):
		makedirs(wopwop_output_path, exist_ok=True)

	if not path.isdir(vtu_output_path):
		makedirs(vtu_output_path, exist_ok=True)

	num_rotors = vehicle.input_state.rotor_inputs.length()
	num_wings = vehicle.input_state.wing_inputs.length()
	print("\n num_wings = ", num_wings, "\n")

	omegas = np.asarray([vehicle.input_state.rotor_inputs[rotor_idx].angular_velocity for rotor_idx in range(num_rotors)])

	log_file.write(f'num_rotors: {num_rotors}\n')
	num_blades = [vehicle.aircraft.rotors[rotor_idx].blades.length() for rotor_idx in range(vehicle.aircraft.rotors.length())]

	# d_psi = computational_parameters['d_psi']

	# dt = d_psi*(math.pi/180.0)/np.max(np.abs(omegas))
	# iter_per_rev = 360/d_psi

	d_psi = computational_parameters['d_psi']

	dt = d_psi*(math.pi/180.0)/np.max(np.abs(omegas))
	iter_per_rev = 360/d_psi

	d_psi = d_psi*np.abs(omegas)/np.max(np.abs(omegas))

	vtk_rotors = [build_base_vtu_rotor(vehicle.aircraft.rotors[rotor_idx]) for rotor_idx in range(num_rotors)]
	vtk_wake = build_base_vtu_wake(vehicle.wake_history.history[0])
	vtk_wing = [build_base_vtu_wing(vehicle.aircraft.wings[w_idx]) for w_idx in range(num_wings)]
	
	#C_T_len = int(round(2.0*math.pi/(dt*max(abs(omegas)))))
	C_T_len = np.round(2.0*math.pi/(dt*np.abs(omegas))).astype(dtype=np.int64)

	num_rotors = vehicle.input_state.rotor_inputs.length()
	average_C_T_arrays = [np.zeros(C_T_len[r_idx]) for r_idx in range(num_rotors)]
	average_F_arrays = np.zeros((3, C_T_len.max()))
	average_C_Mx_arrays = [np.zeros(C_T_len[r_idx]) for r_idx in range(num_rotors)]
	average_C_My_arrays = [np.zeros(C_T_len[r_idx]) for r_idx in range(num_rotors)]
	average_C_Q_arrays = [np.zeros(C_T_len[r_idx]) for r_idx in range(num_rotors)]
	
	average_theta_1c_arrays = [np.zeros(C_T_len[r_idx]) for r_idx in range(num_rotors)]
	average_theta_1s_arrays = [np.zeros(C_T_len[r_idx]) for r_idx in range(num_rotors)]

	average_C_Ts = np.zeros(num_rotors)
	average_C_Mxs = np.zeros(num_rotors)
	average_C_Mys = np.zeros(num_rotors)
	average_C_Qs = np.zeros(num_rotors)
	average_Qs = np.zeros(num_rotors)
	average_Ts = np.zeros(num_rotors)
	average_Fs = np.zeros(3)
	average_Mxs = np.zeros(num_rotors)
	average_Mys = np.zeros(num_rotors)
	average_theta = np.zeros(num_rotors)
	min_theta = np.zeros(num_rotors)
	max_theta = np.zeros(num_rotors)
	theta_rms = np.zeros(num_rotors)
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
	
	trackBWIevents = False
	if "trackBWIevents" in computational_parameters:
		trackBWIevents = computational_parameters["trackBWIevents"]

	log_file.write(f"post_conv_revolutions: {post_conv_revolutions}\n")
	if do_compute:

		bpm_file_header = BPMFileHeader(
			_num_sections = elements,
			_sections_are_uniform = False,
			_includes_section_chord = True,
			_includes_section_length = True,
			_includes_section_te_thickness = True,
			_includes_scetion_te_flow_angle = True,
			_includes_tip_lift_curve_slope = False,
			_includes_section_aoa = True,
			_includes_section_freestream = True
		)
		if (trackBWIevents):
			bwi_file_header = BWIFileHeader(
				_sections_are_uniform = True,
				_includes_section_chord = True,
				_includes_section_length = False,
				_includes_section_freestream = False
			)

		loading_files = [[wopwop_input_files_generator.build_wopwop_loading(rotor, blade, int(round(post_conv_revolutions*iter_per_rev)), naca0012_xsection, wopwop_data_path, acoustics["thickness_noise_flag"]) for blade in rotor.blades] for rotor in vehicle.aircraft.rotors]
		bpm_files = []
		bwi_files = []
		for rotor_idx, rotor in enumerate(vehicle.aircraft.rotors):
			bpm_files.append([])
			bwi_files.append([])
			for blade_idx, blade in enumerate(rotor.blades):
				section_lengths = np.zeros(elements, dtype=np.single)
				r = get_r(rotor.blades[0])
				tallied_length = blade.r_c
				for idx, mid in enumerate(r):
					section_lengths[idx] = 2.0*(mid - tallied_length)
					tallied_length = tallied_length + section_lengths[idx]
					
				section_lengths = rotor.radius*section_lengths

				real_c = rotor.radius*np.asarray(get_chord(blade), dtype=np.single)
				
				bpm_files[rotor_idx].append(wopwop_input_files_generator.build_wopwop_bpm(rotor_idx, blade_idx, int(round(post_conv_revolutions*iter_per_rev)), bpm_file_header, real_c, section_lengths, wopwop_data_path))
				if (trackBWIevents):
					bwi_files[rotor_idx].append(wopwop_input_files_generator.build_wopwop_bwi(rotor_idx, blade_idx, int(round(post_conv_revolutions*iter_per_rev)), bwi_file_header, real_c, section_lengths, wopwop_data_path))
	
	for rotor_idx, rotor in enumerate(vehicle.aircraft.rotors):
		for blade in rotor.blades:
			wopwop_input_files_generator.write_wopwop_geometry(naca0012_xsection, wopwop_data_path, rotor, blade, acoustics["thickness_noise_flag"], omegas[rotor_idx])

	log_file.write("Performing simulation\n")

	trim = False
	moment_trim = False

	trim_algo = 'he'

	trim_mode = TRIM_MODE_COLLECTIVE
	num_trim_groups = num_rotors
	collective_groups = [list(range(num_rotors))]
	collective_force_components = []
	if "trim_mode" in flight_condition:
		if flight_condition["trim_mode"] == "collective":
			trim_mode = TRIM_MODE_COLLECTIVE
		elif flight_condition["trim_mode"] == "shared_collective":
			trim_mode = TRIM_MODE_SHARED_COLLECTIVE
			num_trim_groups = 1
			if "collective_groups" in flight_condition:
				num_trim_groups = len(flight_condition["collective_groups"])
				collective_groups = flight_condition["collective_groups"]

			if "collective_force_components" in flight_condition:
				collective_force_components = flight_condition["collective_force_components"]
			else:
				raise Exception("collective_force_components not specified in parameters file in shared collective trim mode")
			
		elif flight_condition["trim_mode"] == "rpm":
			if "collective_groups" in flight_condition:
				num_trim_groups = len(flight_condition["collective_groups"])
				collective_groups = flight_condition["collective_groups"]
				rotor_signs = [math.copysign(1.0, omegas[collective_groups[trim_group_idx][0]]) for trim_group_idx in range(num_trim_groups)]
			else:
				rotor_signs = [math.copysign(1.0, omegas[rotor_idx]) for rotor_idx in range(num_rotors)]
				collective_groups = [[rotor_idx] for rotor_idx in range(num_rotors)]

			trim_mode = TRIM_MODE_RPM
			needs_trim = [True for _ in range(num_trim_groups)]
			settle_timer = np.zeros(num_rotors)
			omegas = [abs(omega) for omega in omegas]

	radii = [vehicle.aircraft.rotors[rotor_idx].radius for rotor_idx in range(num_rotors)]

	average_theta_arrays = [np.zeros(C_T_len[r_idx]) for r_idx in range(num_trim_groups)]

	if trim_mode == TRIM_MODE_SHARED_COLLECTIVE:
		trim_omegas = [np.asarray([abs(vehicle.input_state.rotor_inputs[r_idx].angular_velocity) for r_idx in collective_groups[trim_group]]).mean() for trim_group in range(num_trim_groups)]
		total_blades = [np.asarray([vehicle.aircraft.rotors[r_idx].blades.length() for r_idx in collective_groups[trim_group]]).sum() for trim_group in range(num_trim_groups)]
		trim_sigmas = [np.asarray([vehicle.aircraft.rotors[r_idx].solidity for r_idx in collective_groups[trim_group]]).mean() for trim_group in range(num_trim_groups)]
		trim_chords = [np.asarray([vehicle.aircraft.rotors[r_idx].blades[0].average_chord for r_idx in collective_groups[trim_group]]).mean() for trim_group in range(num_trim_groups)]
		#trim_radii = [np.asarray([vehicle.aircraft.rotors[r_idx].radius for r_idx in collective_groups[trim_group]]).mean() for trim_group in range(num_trim_groups)]
		trim_radii = [total_blades[trim_group]*trim_chords[trim_group]/trim_sigmas[trim_group] for trim_group in range(num_trim_groups)]
		
		trim_areas = [math.pi*trim_radii[trim_group]**2.0 for trim_group in range(num_trim_groups)]
		print(f'trim_omegas: {trim_omegas}')
		trim_norms = [(atmo.density*trim_areas[trim_group]*trim_radii[trim_group]**2.0*abs(trim_omegas[trim_group])**2.0) for trim_group in range(num_trim_groups)]
		print(f'trim_areas: {trim_areas}')
		print(f'trim_radii: {trim_radii}')
		print(f'trim_norms: {trim_norms}')
	elif trim_mode == TRIM_MODE_RPM:
		trim_norms = [(atmo.density*math.pi*radii[collective_groups[trim_group_idx][0]]**4.0*abs(omegas[collective_groups[trim_group_idx][0]])**2.0) for trim_group_idx in range(num_trim_groups)]

	c_t_bars = None
	c_mx_bars = None
	c_my_bars = None
	if "c_t" in flight_condition:
		trim = True
		c_t_bars = flight_condition["c_t"]

	elif "T" in flight_condition:
		trim = True
		if trim_mode == TRIM_MODE_COLLECTIVE:
			c_t_bars =  [flight_condition["T"][rotor_idx]/(atmo.density*math.pi*vehicle.aircraft.rotors[rotor_idx].radius**4.0*abs(omegas[rotor_idx])**2.0) for rotor_idx in range(num_rotors)]
		elif (trim_mode == TRIM_MODE_SHARED_COLLECTIVE) or (trim_mode == TRIM_MODE_RPM):
			#c_t_bars =  [flight_condition["T"][trim_group] for trim_group in range(num_trim_groups)]
			c_t_bars = [flight_condition["T"][trim_group]/trim_norms[trim_group] for trim_group in range(num_trim_groups)]
	
	

	if "c_mx" in flight_condition and "c_my" in flight_condition:
		moment_trim = True
		c_mx_bars = [flight_condition["c_mx"][rotor_idx]*atmo.density*math.pi*vehicle.aircraft.rotors[rotor_idx].radius**5.0*abs(omegas[rotor_idx])**2.0 for rotor_idx in range(num_rotors)]
		c_my_bars = [flight_condition["c_my"][rotor_idx]*atmo.density*math.pi*vehicle.aircraft.rotors[rotor_idx].radius**5.0*abs(omegas[rotor_idx])**2.0 for rotor_idx in range(num_rotors)]

	elif "Mx" in flight_condition and "My" in flight_condition:
		moment_trim = True
		c_mx_bars = flight_condition["Mx"]
		c_my_bars = flight_condition["My"]

		log_file.write(f'c_mx_bars: {c_mx_bars}, c_my_bars: {c_my_bars}\n')

	if "trim_algo" in flight_condition:
		trim_algo = flight_condition["trim_algo"]
	elif "trim_algo" in computational_parameters:
		trim_algo = computational_parameters["trim_algo"]

	#Ks = 0.4*np.ones(num_trim_groups)
	if trim_mode == TRIM_MODE_RPM:
		Ks = 0.00005*np.ones(num_trim_groups)
		taus = 0.0000001*np.ones(num_trim_groups)
		error = np.ones((num_trim_groups, 3))
		for trim_group_idx in range(num_trim_groups):
			average_C_T_arrays[collective_groups[trim_group_idx][0]][:] = flight_condition["T"][trim_group_idx]/(atmo.density*math.pi*radii[collective_groups[trim_group_idx][0]]**4.0*abs(omegas[collective_groups[trim_group_idx][0]])**2.0)
			#average_Ts[collective_groups[trim_group_idx][0]] = flight_condition["T"][trim_group_idx]
		#Ks = 10*np.ones(num_trim_groups)
		#taus = 20*np.ones(num_trim_groups)
	elif trim_mode == TRIM_MODE_SHARED_COLLECTIVE:
		Ks = 10000*np.ones(num_trim_groups)
		#Ks = 0.01*np.ones(num_trim_groups)
		taus = 0.4*np.ones(num_trim_groups)
	else:
		Ks = 0.4*np.ones(num_trim_groups)
		taus = 4*np.ones(num_trim_groups)
	#taus = 0.4*np.ones(num_trim_groups)
	moment_Ks = np.ones((num_rotors, 2))
	moment_taus = np.ones((num_rotors, 2))

	# HART II aerodas constants
	moment_Ks[0, 0] = 0.3
	moment_Ks[0, 1] = 0.3

	moment_taus[0, 0] = 2.5
	moment_taus[0, 1] = 2.5

	thetas = np.zeros((num_trim_groups, 2))
	betas = np.zeros((num_rotors, max(num_blades), 2))
	
	last_thetas = np.zeros(num_trim_groups)
	moment_thetas = np.zeros((num_rotors, 2*2))

	theta_1s = np.zeros(num_rotors)
	theta_1c = np.zeros(num_rotors)

	rotor_phases = [0 for _ in range(num_rotors)]

	for rotor_idx in range(num_trim_groups):
		last_thetas[rotor_idx] = 10000
		if trim_mode == TRIM_MODE_COLLECTIVE:
			thetas[rotor_idx, 1] = vehicle.input_state.rotor_inputs[rotor_idx].blade_pitches[0]
		elif trim_mode == TRIM_MODE_SHARED_COLLECTIVE:
			thetas[rotor_idx, 1] = vehicle.input_state.rotor_inputs[ collective_groups[rotor_idx][0]].blade_pitches[0]
		if trim_mode == TRIM_MODE_RPM:
			thetas[rotor_idx, 1] = abs(vehicle.input_state.rotor_inputs[rotor_idx].angular_velocity)

	curr_c_ts = np.zeros(num_rotors)
	curr_forces = np.zeros(num_trim_groups)
	#current_ac_forces = Vec4(0)

	curr_c_mxs = np.zeros(num_rotors)
	curr_c_mys = np.zeros(num_rotors)

	last_wake_points = [None for _ in range(num_rotors)]
	last_C_T = [100 for _ in range(num_rotors)]
	last_ac_fz = 10000
	last_C_Mx = [100 for _ in range(num_rotors)]
	last_C_My = [100 for _ in range(num_rotors)]

	wake_l2 = [1000 for _ in range(num_rotors)]
	sim_done = False
	converged = False
	converged_revolutions = 0
	iteration = 0
	acoustic_iteration = 0
	spanwise_element_iteration = 0#[0 for _ in range(num_rotors)]
	aoa_update_iter = 0

	start_time = time.perf_counter_ns()

	theta_3 = 0
	psi_3 = 0

	temp_wake_array_x = np.zeros(wake_lengths[0])
	temp_wake_array_y = np.zeros(wake_lengths[0])
	temp_wake_array_z = np.zeros(wake_lengths[0])
	temp_wake_array_r_c = np.zeros(wake_lengths[0])

	temp_wake_array_x_ts = [[np.zeros(wake_lengths[r_idx]) for _ in range(num_blades[r_idx])] for r_idx in range(num_rotors)]
	temp_wake_array_y_ts = [[np.zeros(wake_lengths[r_idx]) for _ in range(num_blades[r_idx])] for r_idx in range(num_rotors)]
	temp_wake_array_z_ts = [[np.zeros(wake_lengths[r_idx]) for _ in range(num_blades[r_idx])] for r_idx in range(num_rotors)]

	track_wake_element = 'element_trajectories' in results
	track_span_element = 'spanwise_time_series' in results
	track_piv_window = 'piv_window' in results

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

	# wake_idx =  [[[[[[] for _ in range(int(round(post_conv_revolutions*iter_per_rev)) + 1)] for _ in range(max(num_blades))] for _ in range(max(num_blades))] for _ in range(num_rotors)] for _ in range(num_rotors)]
	# bladeSec_idx = [[[[[[] for _ in range(int(round(post_conv_revolutions*iter_per_rev)) + 1)] for _ in range(max(num_blades))] for _ in range(max(num_blades))] for _ in range(num_rotors)] for _ in range(num_rotors)]
	# #wake_miss_dist = [[[[] for _ in range(int(round(post_conv_revolutions*iter_per_rev)) + 1)] for _ in range(max(num_blades))] for _ in range(max(num_blades))]
	# blade_directionVec = [[[[[[] for _ in range(int(round(post_conv_revolutions*iter_per_rev)) + 1)] for _ in range(max(num_blades))] for _ in range(max(num_blades))] for _ in range(num_rotors)] for _ in range(num_rotors)]
	# directionVec_bv = [[[[[[] for _ in range(int(round(post_conv_revolutions*iter_per_rev)) + 1)] for _ in range(max(num_blades))] for _ in range(max(num_blades))] for _ in range(num_rotors)] for _ in range(num_rotors)]
	# vortex_directionVec = [[[[[[] for _ in range(int(round(post_conv_revolutions*iter_per_rev)) + 1)] for _ in range(max(num_blades))] for _ in range(max(num_blades))]for _ in range(num_rotors)] for _ in range(num_rotors)]
	# normalVec = [[[[[[] for _ in range(int(round(post_conv_revolutions*iter_per_rev)) + 1)] for _ in range(max(num_blades))] for _ in range(max(num_blades))] for _ in range(num_rotors)] for _ in range(num_rotors)]
	# missDist = [[[[[[] for _ in range(int(round(post_conv_revolutions*iter_per_rev)) + 1)] for _ in range(max(num_blades))] for _ in range(max(num_blades))] for _ in range(num_rotors)] for _ in range(num_rotors)]
	# gammaSec =[[[[[[] for _ in range(int(round(post_conv_revolutions*iter_per_rev)) + 1)] for _ in range(max(num_blades))] for _ in range(max(num_blades))] for _ in range(num_rotors)] for _ in range(num_rotors)]
	# gammaVortex = [[[[[[] for _ in range(int(round(post_conv_revolutions*iter_per_rev)) + 1)] for _ in range(max(num_blades))] for _ in range(max(num_blades))] for _ in range(num_rotors)] for _ in range(num_rotors)]
	# r_c = [[[[[[] for _ in range(int(round(post_conv_revolutions*iter_per_rev)) + 1)] for _ in range(max(num_blades))] for _ in range(max(num_blades))] for _ in range(num_rotors)] for _ in range(num_rotors)]
	# C_d = [[[[[[] for _ in range(int(round(post_conv_revolutions*iter_per_rev)) + 1)] for _ in range(max(num_blades))] for _ in range(max(num_blades))] for _ in range(num_rotors)] for _ in range(num_rotors)]
	# l = [[[[[[] for _ in range(int(round(post_conv_revolutions*iter_per_rev)) + 1)] for _ in range(max(num_blades))] for _ in range(max(num_blades))] for _ in range(num_rotors)] for _ in range(num_rotors)]
	# secLen =[[[[[[] for _ in range(int(round(post_conv_revolutions*iter_per_rev)) + 1)] for _ in range(max(num_blades))] for _ in range(max(num_blades))] for _ in range(num_rotors)] for _ in range(num_rotors)]
	# #TKE = [[[[] for _ in range(int(round(post_conv_revolutions*iter_per_rev)) + 1)] for _ in range(max(num_blades))] for _ in range(max(num_blades))]
	
	class Interaction:
		def __init__(self):
			self.t = []
			self.blade_idx = []
			self.wake_idx = []
			self.directionVec_bv = []
			self.secIdx = []
			self.vortex_directionVec = []
			self.blade_directionVec = []
			self.normalVec = []
			self.missDist = []
			self.missDist2 = []
			self.gamma = []
			self.Cd = []
			self.r_c = []
			self.l = []
			self.secLen = []
			self.a0 = []
			self.b_e = []
			self.psi = []
			self.ID = []
			self.tipV_idx = []
			self.gamma0 = 0.0
			self.gamma00 = []
			self.Uref = []
			self.L0 = []

	temp_u_p = np.zeros(elements)
	temp_dC_T = np.zeros(elements)
	temp_buffer = np.zeros(elements)
	blade_inflow_distribution = np.zeros((num_rotors, max(num_blades), elements))
	blade_loading_distribution = np.zeros((num_rotors, max(num_blades), elements))

	blade_induced_drag_distribution = np.zeros((num_rotors, max(num_blades), elements))
	blade_profile_drag_distribution = np.zeros((num_rotors, max(num_blades), elements))
	blade_dynamic_aoa_distribution = np.zeros((num_rotors, max(num_blades), elements))

	target_y_slices = []
	piv_slices = []
	piv_window_x = []
	piv_window_z = []
	if track_piv_window:
		piv_slices = results["piv_window"]["y"]
		piv_window_x = results["piv_window"]["x"]
		piv_window_z = results["piv_window"]["z"]
	
	if track_wake_element:
		target_y_slices = results['element_trajectories']

	convergence_rev_multiple = 1
	if "convergence_rev_multiple" in computational_parameters:
		convergence_rev_multiple = computational_parameters["convergence_rev_multiple"]

	span_element_loading = np.zeros((num_rotors, len(target_span_elements), int(convergence_rev_multiple*iter_per_rev)))
	span_element_af_loading = np.zeros((num_rotors, len(target_span_elements), int(convergence_rev_multiple*iter_per_rev)))
	span_element_aoa_eff = np.zeros((num_rotors, len(target_span_elements), int(convergence_rev_multiple*iter_per_rev)))
	span_element_aoa = np.zeros((num_rotors, len(target_span_elements), int(convergence_rev_multiple*iter_per_rev)))
	span_element_up = np.zeros((num_rotors, len(target_span_elements), int(convergence_rev_multiple*iter_per_rev)))
	span_element_inflow_angle = np.zeros((num_rotors, len(target_span_elements), int(convergence_rev_multiple*iter_per_rev)))
	span_element_theta = np.zeros((num_rotors, len(target_span_elements), int(convergence_rev_multiple*iter_per_rev)))
	span_element_gamma = np.zeros((num_rotors, len(target_span_elements), int(convergence_rev_multiple*iter_per_rev)))
	span_element_azimuth = np.zeros((num_rotors, int(convergence_rev_multiple*iter_per_rev)))

	wake_element_index = np.zeros((num_rotors, len(target_y_slices)), dtype=int)
	piv_window_index = np.zeros((num_rotors, max(num_blades), len(piv_slices)), dtype=int)
	wake_element_blade = np.zeros((num_rotors, len(target_y_slices)), dtype=int)
	wake_element_found = [[False for _ in range(len(target_y_slices))] for _ in range(num_rotors)]

	wake_element_trajectory = np.zeros((num_rotors, len(target_y_slices), 2, int(post_conv_revolutions*iter_per_rev) + 1))
	wake_element_core_size = np.zeros((num_rotors, len(target_y_slices), int(post_conv_revolutions*iter_per_rev) + 1))
	
	wake_element_piv = np.zeros((num_rotors, max(num_blades), len(piv_slices), 2, int(post_conv_revolutions*iter_per_rev) + 1))

	actual_wake_history = [wake_lengths[r_idx] if wake_lengths[r_idx]%chunk_size() == 0 else wake_lengths[r_idx] + (chunk_size() - wake_lengths[r_idx]%chunk_size()) for r_idx in range(num_rotors)]
	#wake_trajectory_timehistories = [np.zeros((int((post_conv_revolutions + 1)*iter_per_rev), num_blades[r_idx], 3, actual_wake_history[r_idx])) for r_idx in range(num_rotors)]

	wake_trail_iterations = max(actual_wake_history)

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
	# 	#write_rotor_vtu(f"{vtu_output_path}/rotor", 1000000000, rotor_idx, vtk_rotors[rotor_idx], vehicle.ac_state.rotor_states[rotor_idx], vehicle.input_state.rotor_inputs[rotor_idx], vehicle.aircraft.rotors[rotor_idx])
	# 	write_rotor_vtu(f"{vtu_output_path}/rotor", 1000000000, rotor_idx, vtk_rotors[rotor_idx], vehicle.ac_state.rotor_states[rotor_idx], vehicle.aircraft.rotors[rotor_idx])
	
	convergence_rev_multiple = 1
	if "convergence_rev_multiple" in computational_parameters:
		convergence_rev_multiple = computational_parameters["convergence_rev_multiple"]

	for rotor_idx, rotor in enumerate(vehicle.aircraft.rotors):
		wopwop_motion[rotor.frame.name]["omega"] = vehicle.input_state.rotor_inputs[rotor_idx].angular_velocity

	max_l2 = 1000
	average_step = 0
	if do_compute:

		for m in flight_condition["motion"]:
			if 'blade_element_func' in m:
				if m['blade_element_func'] == 'flapping':
					blade_flapping = lambda az, a=m['cos'], b=m['sin']: flapping_at_azimuth(a, b, 1.0, az)
				elif m['blade_element_func'] == 'pitching':
					elastic_twist = lambda az, a=m['cos'], b=m['sin']: elastic_twist_at_azimuth(a, b, az)

		z_loading = np.zeros(elements, dtype=np.single)
		x_loading = np.zeros(elements, dtype=np.single)
		y_loading = np.zeros(elements, dtype=np.single)
		aoa_array = np.zeros(elements, dtype=np.single)
		u_t = np.zeros(elements, dtype=np.single)
		u_p_array = np.zeros(elements, dtype=np.single)
		u = np.zeros(elements, dtype=np.single)

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

				if (max_l2 is not None) and (max_l2 <= computational_parameters["convergence_criteria"]) and (iteration > (wake_trail_iterations + 2*iter_per_rev)):
					if not converged:
						log_file.write(f"Simulation reached convergence criteria: {max_l2}\n")

					converged = True

			if  convergence_type == 'run_for':
				if iteration/iter_per_rev == flight_condition['run_for']:
					converged = True
					print("converged = ", converged)

			#if iteration % 1 == 0:
			if iteration % iter_per_rev == 0:
				now = time.perf_counter_ns()
				elapsed = now - start_time

				cyclic_inputs = ""
				if moment_trim:
					cyclic_inputs = ''.join([f' θ1c_{rotor_idx}: {average_theta_1c[rotor_idx]*(180.0/math.pi):.4f},' for rotor_idx in range(num_rotors)])+''.join([f' θ1s_{rotor_idx}: {average_theta_1s[rotor_idx]*(180.0/math.pi):.4f},' for rotor_idx in range(num_rotors)])

				if trim_mode == TRIM_MODE_COLLECTIVE:
					log_file.write(
						f'{elapsed/(1000**3):.5f}: rotor rev: {iteration/iter_per_rev:.3f},'
						+f'{average_step/(1000**3):.5f}: '
						+''.join([f' C_T{rotor_idx}: {average_C_Ts[rotor_idx]:.4f},' for rotor_idx in range(num_rotors)])
						+''.join([f' T{rotor_idx}: {average_Ts[rotor_idx]:.4f},' for rotor_idx in range(num_rotors)])
						+''.join([f' P{rotor_idx}: {average_Qs[rotor_idx]:.4f},' for rotor_idx in range(num_rotors)])
						+''.join([f' Mx_{rotor_idx}: {average_Mxs[rotor_idx]:.4f},' for rotor_idx in range(num_rotors)])
						+''.join([f' My_{rotor_idx}: {average_Mys[rotor_idx]:.4f},' for rotor_idx in range(num_rotors)])
						+''.join([f' θ{rotor_idx}: {average_theta[rotor_idx]*(180.0/math.pi):.4f},' for rotor_idx in range(num_rotors)])
						+''.join([f' θmin{rotor_idx}: {min_theta[rotor_idx]*(180.0/math.pi):.4f},' for rotor_idx in range(num_rotors)])
						+''.join([f' θmax{rotor_idx}: {max_theta[rotor_idx]*(180.0/math.pi):.4f},' for rotor_idx in range(num_rotors)])
						+''.join([f' θRMS{rotor_idx}: {theta_rms[rotor_idx]*(180.0/math.pi):.4f},' for rotor_idx in range(num_rotors)])
						+cyclic_inputs
						+''.join([f' χ{rotor_idx}: {average_chis[rotor_idx]:.2f},' for rotor_idx in range(num_rotors)])
						+f' combined_C_T: {average_C_Ts.sum():.5f}, max L_2: {max_l2}\n'
					)
				elif trim_mode == TRIM_MODE_SHARED_COLLECTIVE:
					log_file.write(
						f'{elapsed/(1000**3):.5f}: rotor rev: {iteration/iter_per_rev:.3f},'
						+''.join([f' C_T{rotor_idx}: {average_C_Ts[rotor_idx]:.4f},' for rotor_idx in range(num_rotors)])
						+''.join([f' T{trim_group}: {curr_forces[trim_group]*trim_norms[trim_group]:.4f},' for trim_group in range(num_trim_groups)])
						+''.join([f' P{rotor_idx}: {average_Qs[rotor_idx]:.4f},' for rotor_idx in range(num_rotors)])
						+''.join([f' θ{trim_group}: {average_theta[trim_group]*(180.0/math.pi):.4f},' for trim_group in range(num_trim_groups)])
						+''.join([f' θmin{trim_group}: {min_theta[trim_group]*(180.0/math.pi):.4f},' for trim_group in range(num_trim_groups)])
						+''.join([f' θmax{trim_group}: {max_theta[trim_group]*(180.0/math.pi):.4f},' for trim_group in range(num_trim_groups)])
						+''.join([f' θRMS{trim_group}: {theta_rms[trim_group]*(180.0/math.pi):.4f},' for trim_group in range(num_trim_groups)])
						+cyclic_inputs
						+''.join([f' χ{rotor_idx}: {average_chis[rotor_idx]:.2f},' for rotor_idx in range(num_rotors)])
						+f' max L_2: {max_l2}\n'
						#+f' combined_C_T: {average_C_Ts.sum():.5f}, max L_2: {max_l2}\n'
					)
				elif trim_mode == TRIM_MODE_RPM:
					log_file.write(
						f'{elapsed/(1000**3):.5f}: rotor rev: {iteration/iter_per_rev:.3f},'
						+''.join([f' C_T{rotor_idx}: {average_C_Ts[rotor_idx]:.4f},' for rotor_idx in range(num_rotors)])
						+''.join([f' T{trim_group_idx}: {np.sum([average_Ts[rotor_idx] for rotor_idx in collective_groups[trim_group_idx]]):.4f},' for trim_group_idx in range(num_trim_groups)])
						+''.join([f' P{rotor_idx}: {average_Qs[rotor_idx]:.4f},' for rotor_idx in range(num_rotors)])
						+''.join([f' Ω{rotor_idx}: {vehicle.input_state.rotor_inputs[rotor_idx].angular_velocity*9.5493:.4f},' for rotor_idx in range(num_rotors)])
						+cyclic_inputs
						+''.join([f' χ{rotor_idx}: {average_chis[rotor_idx]:.2f},' for rotor_idx in range(num_rotors)])
						+f' max L_2: {max_l2}\n'
						#+f' combined_C_T: {average_C_Ts.sum():.5f}, max L_2: {max_l2}\n'
					)

				average_step = 0
				# log_file.write(f'theta: {[theta*(180.0/math.pi) for theta in get_theta(vehicle.ac_state.rotor_states[0].blade_states[0])]}\n')
				# log_file.write(f'aoa: {[aoa*(180.0/math.pi) for aoa in get_aoa(vehicle.ac_state.rotor_states[0].blade_states[0])]}\n')
				# log_file.write(f'aoa_eff: {[aoa*(180.0/math.pi) for aoa in get_aoa_eff(vehicle.ac_state.rotor_states[0].blade_states[0])]}\n')
				# log_file.write(f'inflow_angle: {[aoa*(180.0/math.pi) for aoa in get_inflow_angle(vehicle.ac_state.rotor_states[0].blade_states[0])]}\n')
				# log_file.write(f'u_p: {get_u_p(vehicle.ac_state.rotor_states[0].blade_states[0])}\n')
				# log_file.write(f'u_t: {get_u_t(vehicle.ac_state.rotor_states[0].blade_states[0])}\n')
				# log_file.write(f'gamma: {get_gamma(vehicle.ac_state.rotor_states[0].blade_states[0])}\n')
				# log_file.write(f'd_gamma: {get_d_gamma(vehicle.ac_state.rotor_states[0].blade_states[0])}\n')
				# log_file.write(f'dC_T: {get_dC_T(vehicle.ac_state.rotor_states[0].blade_states[0])}\n')

				start_time = now
				log_file.flush()

				if converged and not sim_done:
					if converged_revolutions >= post_conv_revolutions:
						sim_done = True

					converged_revolutions = converged_revolutions + 1

			basic_aircraft_rotor_dynamics(vehicle.input_state, dt)

			for motion_lambda in vehicle.motion_lambdas:
				motion_lambda(vehicle.input_state.rotor_inputs)

			step_start = time.perf_counter_ns()

			# print('before step iteration:', iteration)

			step(vehicle.ac_state, vehicle.aircraft, vehicle.input_state, vehicle.wake_history, atmo, iteration, dt, trackBWIevents, converged)
			average_step = average_step + (time.perf_counter_ns() - step_start)
			# print('after step iteration:', iteration)
			
			# for rotor_idx in range(num_rotors):
			# 	for blade_idx in range(num_rotors):
			# 		for rotor_idx_2 in range(num_rotors):

			# 			len1 = vehicle.wake_history.history[0].rotor_wakes[rotor_idx].interaction_perRotor[rotor_idx_2].blade_vortex_interaction.length()
			# 			len2 = vehicle.wake_history.history[0].rotor_wakes[rotor_idx].interaction_perRotor[rotor_idx_2].blade_vortex_interaction[blade_idx].tip_vortex_interaction.length()
			
			if math.isnan(vehicle.ac_state.rotor_states[0].C_T):
				raise FloatingPointError("Simulation failed: Encountered NaN in rotor thrust")
			
			average_F_arrays[0,iteration % C_T_len.max()] = vehicle.ac_state.forces[0]
			average_F_arrays[1,iteration % C_T_len.max()] = vehicle.ac_state.forces[1]
			average_F_arrays[2,iteration % C_T_len.max()] = vehicle.ac_state.forces[2]

			average_Fs[0] = average_F_arrays[0,:].mean()
			average_Fs[1] = average_F_arrays[1,:].mean()
			average_Fs[2] = average_F_arrays[2,:].mean()

			for rotor_idx, rotor_state in enumerate(vehicle.ac_state.rotor_states):
				
				average_C_T_arrays[rotor_idx][iteration % C_T_len[rotor_idx]] = rotor_state.C_T
				average_C_Q_arrays[rotor_idx][iteration % C_T_len[rotor_idx]] = rotor_state.C_Q
				average_C_Mx_arrays[rotor_idx][iteration % C_T_len[rotor_idx]] = rotor_state.C_Mx
				average_C_My_arrays[rotor_idx][iteration % C_T_len[rotor_idx]] = rotor_state.C_My
				average_chi_arrays[rotor_idx][iteration % C_T_len[rotor_idx]] = vehicle.inflows[rotor_idx].wake_skew()*180.0/math.pi
				
				if trim_mode == TRIM_MODE_COLLECTIVE:
					average_theta_arrays[rotor_idx][iteration % C_T_len[rotor_idx]] = thetas[rotor_idx,1]
				average_theta_1c_arrays[rotor_idx][iteration % C_T_len[rotor_idx]] = theta_1c[rotor_idx]
				average_theta_1s_arrays[rotor_idx][iteration % C_T_len[rotor_idx]] = theta_1s[rotor_idx]

				average_C_Ts[rotor_idx] = np.sum(average_C_T_arrays[rotor_idx])/C_T_len[rotor_idx]
				average_C_Mxs[rotor_idx] = np.sum(average_C_Mx_arrays[rotor_idx])/C_T_len[rotor_idx]
				average_C_Mys[rotor_idx] = np.sum(average_C_My_arrays[rotor_idx])/C_T_len[rotor_idx]
				average_C_Qs[rotor_idx] = np.sum(average_C_Q_arrays[rotor_idx])/C_T_len[rotor_idx]

				if trim_mode == TRIM_MODE_COLLECTIVE:
					average_theta[rotor_idx] = np.sum(average_theta_arrays[rotor_idx])/C_T_len[rotor_idx]
					min_theta[rotor_idx] = average_theta_arrays[rotor_idx].min()
					max_theta[rotor_idx] = average_theta_arrays[rotor_idx].max()
					theta_rms[rotor_idx] = np.sqrt(np.power(average_theta_arrays[rotor_idx], 2.0).mean())

				average_theta_1c[rotor_idx] = np.sum(average_theta_1c_arrays[rotor_idx])/C_T_len[rotor_idx]
				average_theta_1s[rotor_idx] = np.sum(average_theta_1s_arrays[rotor_idx])/C_T_len[rotor_idx]

				average_Ts[rotor_idx] = average_C_Ts[rotor_idx]*atmo.density*math.pi*radii[rotor_idx]**4.0*abs(omegas[rotor_idx])**2.0
				average_Mxs[rotor_idx] = average_C_Mxs[rotor_idx]*atmo.density*math.pi*radii[rotor_idx]**5.0*abs(omegas[rotor_idx])**2.0
				average_Mys[rotor_idx] = -average_C_Mys[rotor_idx]*atmo.density*math.pi*radii[rotor_idx]**5.0*abs(omegas[rotor_idx])**2.0
				average_Qs[rotor_idx] = average_C_Qs[rotor_idx]*atmo.density*math.pi*radii[rotor_idx]**5.0*abs(omegas[rotor_idx])**3.0

				average_chis[rotor_idx] = np.sum(average_chi_arrays[rotor_idx])/C_T_len[rotor_idx]

			if trim_mode == TRIM_MODE_SHARED_COLLECTIVE:
				for trim_group in range(num_trim_groups):
					average_theta_arrays[trim_group][iteration % C_T_len[trim_group]] = thetas[trim_group,1]
					average_theta[trim_group] = np.sum(average_theta_arrays[trim_group])/C_T_len[trim_group]
					min_theta[trim_group] = average_theta_arrays[trim_group].min()
					max_theta[trim_group] = average_theta_arrays[trim_group].max()
					theta_rms[trim_group] = np.sqrt(np.power(average_theta_arrays[trim_group], 2.0).mean())

			if moment_trim:
				for rotor_idx, rotor_state in enumerate(vehicle.ac_state.rotor_states):
					curr_c_mxs[rotor_idx] = rotor_state.C_Mx*atmo.density*math.pi*radii[rotor_idx]**5.0*abs(omegas[rotor_idx])**2.0
					curr_c_mys[rotor_idx] = rotor_state.C_My*atmo.density*math.pi*radii[rotor_idx]**5.0*abs(omegas[rotor_idx])**2.0
				
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
					N = atmo.density*math.pi*radii[rotor_idx]**5.0*abs(omegas[rotor_idx])**2.0

					moment_thetas[rotor_idx,:] = moment_thetas[rotor_idx,:] + dt*abs(omegas[rotor_idx])*moment_trim_ode(moment_thetas[rotor_idx,:], moment_taus[rotor_idx,:], c_mx_bars[rotor_idx], c_my_bars[rotor_idx], moment_Ks[rotor_idx,:], curr_c_mxs[rotor_idx], -curr_c_mys[rotor_idx], N, rotor_idx)
					theta_1s[rotor_idx] = moment_thetas[rotor_idx,1]
					theta_1c[rotor_idx] = moment_thetas[rotor_idx,3]

			if trim_mode == TRIM_MODE_SHARED_COLLECTIVE:
				for trim_group in range(num_trim_groups):
					curr_forces[trim_group] = average_Fs[collective_force_components[trim_group]]/trim_norms[trim_group]
					
			if trim:
				for trim_group_idx in range(num_trim_groups):
					if trim_mode == TRIM_MODE_RPM:
						c_t_bars[trim_group_idx] = flight_condition["T"][trim_group_idx]#/(atmo.density*math.pi*radii[trim_group_idx]**4.0*abs(omegas[trim_group_idx])**2.0)
						#c_t_bars[trim_group_idx] = flight_condition["T"][trim_group_idx]/(atmo.density*math.pi*radii[collective_groups[trim_group_idx][0]]**4.0*abs(omegas[collective_groups[trim_group_idx][0]])**2.0)
						v = [vehicle.ac_state.rotor_states[rotor_idx].C_T*(atmo.density*math.pi*radii[rotor_idx]**4.0*abs(omegas[rotor_idx])**2.0) for rotor_idx in collective_groups[trim_group_idx]]
						#v = [vehicle.ac_state.rotor_states[rotor_idx].C_T for rotor_idx in collective_groups[trim_group_idx]]
						#print(f"trim_group_idx: {trim_group_idx} collective_groups[trim_group_idx]: {collective_groups[trim_group_idx]} v: {v}")
						curr_c_ts[trim_group_idx] = np.sum(v)
						#curr_c_ts[trim_group_idx] = np.sum([vehicle.ac_state.rotor_states[rotor_idx].C_T for rotor_idx in collective_groups[trim_group_idx]])
					else:
						curr_c_ts[trim_group_idx] = vehicle.ac_state.rotor_states[trim_group_idx].C_T

				def trim_ode(theta, tau, c_t_bar, K, curr_c_t):
					return np.asarray([
						1/tau*(K*6/(2.0*math.pi*vehicle.aircraft.rotors[0].solidity)*(c_t_bar - curr_c_t) - theta[0]),
						theta[0]
					])

				if trim_mode == TRIM_MODE_COLLECTIVE:
					for rotor_idx, rotor_state in enumerate(vehicle.ac_state.rotor_states):
						if trim_algo == 'he':
							thetas[rotor_idx,:] = thetas[rotor_idx,:] + dt*abs(omegas[rotor_idx])*trim_ode(thetas[rotor_idx,:], taus[rotor_idx], c_t_bars[rotor_idx], Ks[rotor_idx], curr_c_ts[rotor_idx])
						elif trim_algo == 'lympany':
							thetas[rotor_idx,1] = thetas[rotor_idx,1] + 3.0/(math.pi*vehicle.aircraft.rotors[rotor_idx].solidity)*(c_t_bars[rotor_idx] - curr_c_ts[rotor_idx])
				elif trim_mode == TRIM_MODE_SHARED_COLLECTIVE:
					for trim_group in range(num_trim_groups):
						if trim_algo == 'he':
							print(f'c_t_bars[{trim_group}]: {c_t_bars[trim_group]}, curr_forces[{trim_group}]: {curr_forces[trim_group]} delta: {c_t_bars[trim_group] - curr_forces[trim_group]} thetas[{trim_group}, 1]: {thetas[trim_group, 1]*(180.0/math.pi)}')
						elif trim_algo == 'lympany':
							thetas[trim_group, 1] = thetas[trim_group, 1] + 3.0/(math.pi*vehicle.aircraft.rotors[collective_groups[trim_group][0]].solidity)*(c_t_bars[trim_group] - curr_forces[trim_group])
							print(f'c_t_bars[{trim_group}]: {c_t_bars[trim_group]}, curr_forces[{trim_group}]: {curr_forces[trim_group]} delta: {c_t_bars[trim_group] - curr_forces[trim_group]} thetas[{trim_group}, 1]: {thetas[trim_group, 1]*(180.0/math.pi)}')
				elif trim_mode == TRIM_MODE_RPM:
					if (iteration % 1 == 0) and not converged:
						for trim_group_idx in range(num_trim_groups):
							if needs_trim[trim_group_idx]:
								error[trim_group_idx, 2] = error[trim_group_idx, 1]
								error[trim_group_idx, 1] = error[trim_group_idx, 0]
								error[trim_group_idx, 0] = c_t_bars[trim_group_idx] - curr_c_ts[trim_group_idx]

								de = 1/(2*dt)*(3.0*error[trim_group_idx, 0] - 4*error[trim_group_idx, 1] + error[trim_group_idx, 2])

								thetas[trim_group_idx,1] = thetas[trim_group_idx,1] + Ks[trim_group_idx]*error[trim_group_idx, 0] + de*taus[trim_group_idx]

			if trim_mode == TRIM_MODE_COLLECTIVE:
				for rotor_idx, rotor_state in enumerate(vehicle.ac_state.rotor_states):
					# print(f'vehicle.ac_state.rotor_states.length(): {vehicle.ac_state.rotor_states.length()}')
					# print(f'rotor {rotor_idx}')
					# print(f'num_blades: {num_blades[rotor_idx]}')
					for blade_idx in range(num_blades[rotor_idx]):
						blade_azimuth = vehicle.input_state.rotor_inputs[rotor_idx].azimuth + vehicle.aircraft.rotors[rotor_idx].blades[blade_idx].azimuth_offset
						cos_azimuth = math.cos(blade_azimuth)
						sin_azimuth = math.sin(blade_azimuth)

						sin3_azimuth = math.cos(3.0*(blade_azimuth) - (psi_3 - math.pi))
						#sin3_azimuth = math.cos(3.0*(blade_azimuth) - (psi_3))
						#sin3_azimuth = math.sin(3.0*(blade_azimuth) - psi_3)

						vehicle.input_state.rotor_inputs[rotor_idx].blade_pitches[blade_idx] = thetas[rotor_idx, 1] + theta_1c[rotor_idx]*cos_azimuth + theta_1s[rotor_idx]*sin_azimuth + theta_3*sin3_azimuth
						#print(f'trimming rotor {rotor_idx} blade {blade_idx} to {vehicle.input_state.rotor_inputs[rotor_idx].blade_pitches[blade_idx]}')
						if elastic_twist and (rotor_idx == 0):
							#vehicle.input_state.rotor_inputs[rotor_idx].blade_pitches[blade_idx] = thetas[rotor_idx, 1] + theta_1c[rotor_idx]*cos_azimuth + theta_1s[rotor_idx]*sin_azimuth + theta_3*sin3_azimuth
							#vehicle.input_state.rotor_inputs[rotor_idx].blade_pitches[blade_idx] = vehicle.input_state.rotor_inputs[rotor_idx].blade_pitches[blade_idx] + elastic_twist(blade_azimuth)
							vehicle.input_state.rotor_inputs[rotor_idx].blade_pitches[blade_idx] = vehicle.input_state.rotor_inputs[rotor_idx].blade_pitches[blade_idx] + elastic_twist(blade_azimuth - math.pi)
						# else:
						# 	#vehicle.input_state.rotor_inputs[rotor_idx].blade_pitches[blade_idx] = math.copysign(1, omegas[rotor_idx])*(thetas[rotor_idx, 1] + theta_1c[rotor_idx]*cos_azimuth + theta_1s[rotor_idx]*sin_azimuth + theta_3*sin3_azimuth)
						# 	vehicle.input_state.rotor_inputs[rotor_idx].blade_pitches[blade_idx] = thetas[rotor_idx, 1] + theta_1c[rotor_idx]*cos_azimuth + theta_1s[rotor_idx]*sin_azimuth + theta_3*sin3_azimuth
			elif trim_mode == TRIM_MODE_SHARED_COLLECTIVE:
				for trim_group in range(num_trim_groups):
					for rotor_idx in collective_groups[trim_group]:
						for blade_idx in range(num_blades[rotor_idx]):
							set_blade_pitch(vehicle.input_state, rotor_idx, blade_idx, thetas[trim_group, 1])

			elif trim_mode == TRIM_MODE_RPM:
				if (iteration % 1 == 0) and not converged:
					for trim_group_idx in range(num_trim_groups):
						if needs_trim[trim_group_idx]:

							for rotor_idx in collective_groups[trim_group_idx]:
								vehicle.input_state.rotor_inputs[rotor_idx].angular_velocity = rotor_signs[trim_group_idx]*thetas[trim_group_idx, 1]
								omegas[rotor_idx] = thetas[trim_group_idx, 1]

							if iteration > wake_trail_iterations:
								c_t_delta = np.sum([average_Ts[rotor_idx] for rotor_idx in collective_groups[trim_group_idx]]) - c_t_bars[trim_group_idx]

								if abs(c_t_delta/c_t_bars[trim_group_idx]) <= 0.05:

									vehicle.input_state.rotor_inputs[rotor_idx].angular_accel = 0
									settle_timer[trim_group_idx] = 0
						else:
							settle_timer[trim_group_idx] = settle_timer[trim_group_idx] + 1

							c_t_delta = np.sum([average_Ts[rotor_idx] for rotor_idx in collective_groups[trim_group_idx]]) - c_t_bars[trim_group_idx]
							if (abs(c_t_delta ) > 5.0e-5) and (settle_timer[trim_group_idx] >= 2*iter_per_rev):
								print(f"Starting trim for rotor {trim_group_idx}")
								needs_trim[trim_group_idx] = True

				#if not converged:
				dt = computational_parameters['d_psi']*(math.pi/180.0)/np.max(np.abs(np.asarray(omegas)))
			numBladesArray = np.zeros(num_rotors, dtype=int)
			for rotor_idx, rotor_state in enumerate(vehicle.ac_state.rotor_states):
				
				numBladesArray[rotor_idx] = rotor_state.blade_states.length()
				# print(numBladesArray)
				if ((rotor_idx == 1) and (vehicle.name == "helinovi")) or ((rotor_idx == 0) and (vehicle.name == "helinovi_tr")):
					net_blade_moment = rotor_state.blade_states[0].C_My
					# for _, blade_state in enumerate(rotor_state.blade_states):
					# 	net_blade_moment = net_blade_moment + blade_state.C_My

					gamma = flight_condition["lock_number"][rotor_idx]
					delta_3 = flight_condition["delta_3"]*(math.pi/180.0)

					def flapping_ode(beta):
						return np.asarray([
							beta[1],
							gamma*net_blade_moment - gamma/8.0*beta[1] - (1 + gamma/8.0*math.tan(delta_3))*beta[0]
						])

					betas[rotor_idx,0,:] = betas[rotor_idx,0,:] + dt*abs(omegas[rotor_idx])*flapping_ode(betas[rotor_idx,0,:])

					for blade_idx in range(rotor_state.blade_states.length()):
						if blade_idx == 0:
							vehicle.aircraft.rotors[rotor_idx].blades[blade_idx].frame.parent.set_rotation(Vec3([math.sin(delta_3), math.cos(delta_3), 0]), -betas[rotor_idx, 0, 0])
							vehicle.input_state.rotor_inputs[rotor_idx].blade_flapping_rate[blade_idx] = -betas[rotor_idx, 0, 1]
							vehicle.input_state.rotor_inputs[rotor_idx].blade_pitches[blade_idx] = vehicle.input_state.rotor_inputs[rotor_idx].blade_pitches[blade_idx] + math.tan(delta_3)*betas[rotor_idx, 0, 0]
						else:
							vehicle.aircraft.rotors[rotor_idx].blades[blade_idx].frame.parent.set_rotation(Vec3([math.sin(delta_3), math.cos(delta_3), 0]), betas[rotor_idx, 0, 0])
							vehicle.input_state.rotor_inputs[rotor_idx].blade_flapping_rate[blade_idx] = betas[rotor_idx, 0, 1]
							vehicle.input_state.rotor_inputs[rotor_idx].blade_pitches[blade_idx] = vehicle.input_state.rotor_inputs[rotor_idx].blade_pitches[blade_idx] - math.tan(delta_3)*betas[rotor_idx, 0, 0]

				elif blade_flapping is not None:
					for blade_idx in range(num_blades[rotor_idx]):
						blade_azimuth = vehicle.input_state.rotor_inputs[rotor_idx].azimuth + vehicle.aircraft.rotors[rotor_idx].blades[blade_idx].azimuth_offset
						(h, h_star) = blade_flapping(blade_azimuth - math.pi)
						#(h, h_star) = blade_flapping(blade_azimuth)
						vehicle.input_state.rotor_inputs[rotor_idx].blade_flapping_rate[blade_idx] = h_star

			loading_data.time = dt*acoustic_iteration

			if converged:
				# print('converged: acoustic_iteration:', acoustic_iteration);
				if write_wake and (converged_revolutions >= (post_conv_revolutions - 1)):
					for rotor_idx in range(num_rotors):
						print("writing rotor and wake vtu")
						write_rotor_vtu(f"{vtu_output_path}/rotor", acoustic_iteration, rotor_idx, vtk_rotors[rotor_idx], vehicle.ac_state.rotor_states[rotor_idx], vehicle.aircraft.rotors[rotor_idx])
						write_wake_vtu(f"{vtu_output_path}/wake", acoustic_iteration, vtk_wake, vehicle.wake_history.history[0])
					
					for wing_idx in range(num_wings):
						print("writing wing vtu")
						write_wing_vtu(f"{vtu_output_path}/wing", acoustic_iteration, wing_idx, vtk_wing[wing_idx], vehicle.ac_state.wing_states[wing_idx], vehicle.aircraft.wings[wing_idx])

				if (spanwise_element_iteration >= convergence_rev_multiple*iter_per_rev) and start_recording:
					done_recording = True

				if iteration%(convergence_rev_multiple*iter_per_rev) == 0 and not done_recording:
					start_recording = True

				if acoustic_iteration == 0:
					for rotor_idx in range(num_rotors):
						rotor_phases[rotor_idx] = vehicle.input_state.rotor_inputs[rotor_idx].azimuth
				

				for rotor_idx, rotor in enumerate(vehicle.ac_state.rotor_states):

					wake_idx =  [[[[[] for _ in range(max(num_blades))] for _ in range(max(num_blades))] for _ in range(num_rotors)] for _ in range(num_rotors)]
					bladeSec_idx = [[[[[] for _ in range(max(num_blades))] for _ in range(max(num_blades))] for _ in range(num_rotors)] for _ in range(num_rotors)]
					blade_directionVec = [[[[[] for _ in range(max(num_blades))] for _ in range(max(num_blades))] for _ in range(num_rotors)] for _ in range(num_rotors)]
					directionVec_bv = [[[[[] for _ in range(max(num_blades))] for _ in range(max(num_blades))] for _ in range(num_rotors)] for _ in range(num_rotors)]
					vortex_directionVec = [[[[[] for _ in range(max(num_blades))] for _ in range(max(num_blades))]for _ in range(num_rotors)] for _ in range(num_rotors)]
					normalVec = [[[[[] for _ in range(max(num_blades))] for _ in range(max(num_blades))] for _ in range(num_rotors)] for _ in range(num_rotors)]
					missDist = [[[[[] for _ in range(max(num_blades))] for _ in range(max(num_blades))] for _ in range(num_rotors)] for _ in range(num_rotors)]
					gammaSec =[[[[[] for _ in range(max(num_blades))] for _ in range(max(num_blades))] for _ in range(num_rotors)] for _ in range(num_rotors)]
					gammaVortex = [[[[[] for _ in range(max(num_blades))] for _ in range(max(num_blades))] for _ in range(num_rotors)] for _ in range(num_rotors)]
					r_c = [[[[[] for _ in range(max(num_blades))] for _ in range(max(num_blades))] for _ in range(num_rotors)] for _ in range(num_rotors)]
					C_d = [[[[[] for _ in range(max(num_blades))] for _ in range(max(num_blades))] for _ in range(num_rotors)] for _ in range(num_rotors)]
					l = [[[[[] for _ in range(max(num_blades))] for _ in range(max(num_blades))] for _ in range(num_rotors)] for _ in range(num_rotors)]
					secLen =[[[[[] for _ in range(max(num_blades))] for _ in range(max(num_blades))] for _ in range(num_rotors)] for _ in range(num_rotors)]

					blade_idx = 2
					if (rotor_idx > 0) or ((rotor_idx == 0) and (vehicle.name == "helinovi_tr")):
						blade_idx = 0

					if start_recording and not done_recording:
						for t_idx in range(len(target_span_elements)):
							dC_l = vehicle.ac_state.rotor_states[rotor_idx].blade_states[blade_idx].chunks[target_span_chunk_index[rotor_idx][t_idx]].dC_l[target_span_element_index[rotor_idx][t_idx]]
							dC_L = vehicle.ac_state.rotor_states[rotor_idx].blade_states[blade_idx].chunks[target_span_chunk_index[rotor_idx][t_idx]].dC_L[target_span_element_index[rotor_idx][t_idx]]
							aoa_eff = vehicle.ac_state.rotor_states[rotor_idx].blade_states[blade_idx].chunks[target_span_chunk_index[rotor_idx][t_idx]].aoa_eff[target_span_element_index[rotor_idx][t_idx]]
							aoa = vehicle.ac_state.rotor_states[rotor_idx].blade_states[blade_idx].chunks[target_span_chunk_index[rotor_idx][t_idx]].aoa[target_span_element_index[rotor_idx][t_idx]]
							u_p = vehicle.ac_state.rotor_states[rotor_idx].blade_states[blade_idx].chunks[target_span_chunk_index[rotor_idx][t_idx]].u_p[target_span_element_index[rotor_idx][t_idx]]
							inflow_angle = vehicle.ac_state.rotor_states[rotor_idx].blade_states[blade_idx].chunks[target_span_chunk_index[rotor_idx][t_idx]].inflow_angle[target_span_element_index[rotor_idx][t_idx]]
							theta = vehicle.ac_state.rotor_states[rotor_idx].blade_states[blade_idx].chunks[target_span_chunk_index[rotor_idx][t_idx]].theta[target_span_element_index[rotor_idx][t_idx]]
							gamma = vehicle.ac_state.rotor_states[rotor_idx].blade_states[blade_idx].chunks[target_span_chunk_index[rotor_idx][t_idx]].gamma[target_span_element_index[rotor_idx][t_idx]]

							span_element_af_loading[rotor_idx, t_idx, spanwise_element_iteration] = dC_l
							span_element_loading[rotor_idx, t_idx, spanwise_element_iteration] = dC_L*atmo.density*math.pi*radii[rotor_idx]**3.0*abs(omegas[rotor_idx])**2.0/(0.5*atmo.density*flight_condition["sos"]**2.0*vehicle.aircraft.rotors[rotor_idx].blades[blade_idx].average_chord)
							span_element_aoa_eff[rotor_idx, t_idx, spanwise_element_iteration] = aoa_eff*(180.0/math.pi)
							span_element_aoa[rotor_idx, t_idx, spanwise_element_iteration] = aoa*(180.0/math.pi)
							span_element_up[rotor_idx, t_idx, spanwise_element_iteration] = u_p
							span_element_inflow_angle[rotor_idx, t_idx, spanwise_element_iteration] = inflow_angle*(180.0/math.pi)
							span_element_theta[rotor_idx, t_idx, spanwise_element_iteration] = theta*(180.0/math.pi)
							span_element_gamma[rotor_idx, t_idx, spanwise_element_iteration] = gamma
							span_element_azimuth[rotor_idx, spanwise_element_iteration] = spanwise_element_iteration*d_psi[rotor_idx]
						#spanwise_element_iteration[rotor_idx] = spanwise_element_iteration[rotor_idx] + 1

					interaction = [Interaction() for _ in range(max(num_blades))]
					perpInteraction = [Interaction() for _ in range(max(num_blades))]
					perpInteraction_perBlade = [Interaction() for _ in range(max(num_blades))]

					for blade_idx, blade in enumerate(rotor.blade_states):

						#blade_azimuth = vehicle.input_state.rotor_inputs[rotor_idx].azimuth + vehicle.aircraft.rotors[rotor_idx].blades[blade_idx].azimuth_offset
						blade_azimuth = get_blade_azimuth(vehicle.aircraft, vehicle.input_state, rotor_idx, blade_idx)
						
						cos_azimuth = math.cos(blade_azimuth)
						sin_azimuth = math.sin(blade_azimuth)

						sin3_azimuth = math.cos(3.0*(blade_azimuth) - (psi_3 - math.pi))

						collective_pitch_array[rotor_idx, acoustic_iteration] = thetas[rotor_idx, 1]#get_blade_pitch(vehicle.input_state, rotor_idx, blade_idx)
						#collective_pitch_array[rotor_idx, acoustic_iteration] = vehicle.input_state.rotor_inputs[rotor_idx].blade_pitches[blade_idx] # thetas[rotor_idx, 1]
						sin_pitch_array[rotor_idx, blade_idx, acoustic_iteration] = theta_1s[rotor_idx]
						cos_pitch_array[rotor_idx, blade_idx, acoustic_iteration] = theta_1c[rotor_idx]
						hhc_pitch_array[rotor_idx, blade_idx, acoustic_iteration] = theta_3*sin3_azimuth

						blade_twist_azimuth[rotor_idx, blade_idx, acoustic_iteration] = acoustic_iteration*d_psi[rotor_idx]
						if elastic_twist:
							elastic_twist_array[rotor_idx, blade_idx, acoustic_iteration] = elastic_twist(blade_azimuth - math.pi)

						blade_twist_array[rotor_idx, blade_idx, acoustic_iteration] = get_blade_pitch(vehicle.input_state, rotor_idx, blade_idx)
						#blade_twist_array[rotor_idx, blade_idx, acoustic_iteration] = vehicle.input_state.rotor_inputs[rotor_idx].blade_pitches[blade_idx]

						#if (blade_azimuth < (0.5*math.pi + 0.75*(math.pi/180.0))) and (blade_azimuth > (0.5*math.pi - 0.75*(math.pi/180.0))):
						if (blade_azimuth < (1.5*math.pi + 0.75*(math.pi/180.0))) and (blade_azimuth > (1.5*math.pi - 0.75*(math.pi/180.0))):
							fill_dynamic_u_pd(vehicle.ac_state.rotor_states[rotor_idx].blade_states[blade_idx], temp_buffer)
							blade_inflow_distribution[rotor_idx, blade_idx, :] = temp_buffer

							fill_dC_Td(vehicle.ac_state.rotor_states[rotor_idx].blade_states[blade_idx], temp_buffer)
							blade_loading_distribution[rotor_idx, blade_idx, :] = temp_buffer

							fill_dynamic_dC_Db_profile(vehicle.ac_state.rotor_states[rotor_idx].blade_states[blade_idx], temp_buffer)
							blade_profile_drag_distribution[rotor_idx, blade_idx, :] = temp_buffer

							fill_dynamic_dC_Db_induced(vehicle.ac_state.rotor_states[rotor_idx].blade_states[blade_idx], temp_buffer)
							blade_induced_drag_distribution[rotor_idx, blade_idx, :] = temp_buffer

							fill_aoa_effd(vehicle.ac_state.rotor_states[rotor_idx].blade_states[blade_idx], temp_buffer)
							blade_dynamic_aoa_distribution[rotor_idx, blade_idx, :] = temp_buffer

						if blade_flapping is not None:
							(h, h_star) = blade_flapping(blade_azimuth - math.pi)
							blade_flapping_array[rotor_idx, blade_idx, acoustic_iteration] = h
							blade_flapping_der_array[rotor_idx, blade_idx, acoustic_iteration] = h_star

						if track_wake_element:
							for t_idx, target_y_slice in enumerate(target_y_slices):
								fill_wake_y_component(vehicle.wake_history.history[0].rotor_wakes[rotor_idx].tip_vortices[blade_idx], temp_wake_array_y)

								if (abs(temp_wake_array_y[0] - (radii[rotor_idx]*target_y_slice)) < 0.05) and (not wake_element_found[rotor_idx][t_idx]) and (((blade.azimuth <= 1.02*0.5*math.pi) and (blade.azimuth >= -1.02*0.5*math.pi)) or ((blade.azimuth <= 1.02*5.0/2.0*math.pi) and (blade.azimuth >= 1.02*3.0/2.0*math.pi))):
									num_chunks = blade.chunks.len()
									log_file.write(f"Found wake element to track for slice {target_y_slice} at position: {blade.chunks[num_chunks - 1].x[7]}, {temp_wake_array_y[0]}, {blade.chunks[num_chunks - 1].z[7]}\n")
									wake_element_found[rotor_idx][t_idx] = True
									wake_element_blade[rotor_idx, t_idx] = blade_idx

						z_loading = z_loading.astype(np.float32)

						fill_dC_Nf(blade, z_loading)
						fill_dC_cf(blade, x_loading)

						z_loading = -z_loading*atmo.density*math.pi*radii[rotor_idx]**3.0*abs(omegas[rotor_idx])**2.0
						#x_loading = -x_loading*atmo.density*math.pi*radii[rotor_idx]**3.0*abs(omegas[rotor_idx])**2.0

						loading_data.set_z_loading_array(z_loading.astype(np.float32))
						loading_data.set_y_loading_array(y_loading)
						loading_data.set_x_loading_array(x_loading)
						append_loading_data(loading_files[rotor_idx][blade_idx], loading_data)

						real_c = radii[rotor_idx]*np.asarray(get_chord(vehicle.aircraft.rotors[rotor_idx].blades[blade_idx]), dtype=np.single)
						# Nitya, 09.16
						
						if(trackBWIevents):
							# print('1. acoustic_iteration:', acoustic_iteration)
							for rotor_idx_2 in range(num_rotors):
								#print(vehicle.wake_history.history[0].rotor_wakes[rotor_idx].interaction_perRotor[rotor_idx_2].blade_vortex_interaction[blade_idx].tip_vortex_interaction.length())
								for blade_idx2 in range(0,vehicle.wake_history.history[0].rotor_wakes[rotor_idx_2].interaction_perRotor[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction.length()):
									# wake_idx[rotor_idx][rotor_idx_2][blade_idx][blade_idx2].append( get_interactionPt_wake_idx(vehicle.wake_history.history[0].rotor_wakes[rotor_idx_2].interaction_perRotor[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction[blade_idx2]))
									# bladeSec_idx[rotor_idx][rotor_idx_2][blade_idx][blade_idx2].append(get_interactionPt_bladeSec_idx(vehicle.wake_history.history[0].rotor_wakes[rotor_idx_2].interaction_perRotor[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction[blade_idx2]))
									# blade_directionVec[rotor_idx][rotor_idx_2][blade_idx][blade_idx2].append(get_interaction_point_r_blade(vehicle.wake_history.history[0].rotor_wakes[rotor_idx_2].interaction_perRotor[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction[blade_idx2]))
									# directionVec_bv[rotor_idx][rotor_idx_2][blade_idx][blade_idx2].append(get_interaction_point_r_blade_v(vehicle.wake_history.history[0].rotor_wakes[rotor_idx_2].interaction_perRotor[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction[blade_idx2]))
									# vortex_directionVec[rotor_idx][rotor_idx_2][blade_idx][blade_idx2].append(get_interaction_point_r_vortex(vehicle.wake_history.history[0].rotor_wakes[rotor_idx_2].interaction_perRotor[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction[blade_idx2]))
									# normalVec[rotor_idx][rotor_idx_2][blade_idx][blade_idx2].append(get_interaction_point_blade_normal(vehicle.wake_history.history[0].rotor_wakes[rotor_idx_2].interaction_perRotor[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction[blade_idx2]))
									# missDist[rotor_idx][rotor_idx_2][blade_idx][blade_idx2].append(get_interaction_point_miss_dist(vehicle.wake_history.history[0].rotor_wakes[rotor_idx_2].interaction_perRotor[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction[blade_idx2]))
									# gammaSec[rotor_idx][rotor_idx_2][blade_idx][blade_idx2].append(get_interaction_point_gamma_sec(vehicle.wake_history.history[0].rotor_wakes[rotor_idx_2].interaction_perRotor[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction[blade_idx2]))
									# gammaVortex[rotor_idx][rotor_idx_2][blade_idx][blade_idx2].append(get_interaction_point_gamma_w(vehicle.wake_history.history[0].rotor_wakes[rotor_idx_2].interaction_perRotor[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction[blade_idx2]))
									# r_c[rotor_idx][rotor_idx_2][blade_idx][blade_idx2].append(get_interaction_point_r_c(vehicle.wake_history.history[0].rotor_wakes[rotor_idx_2].interaction_perRotor[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction[blade_idx2]))
									# C_d[rotor_idx][rotor_idx_2][blade_idx][blade_idx2].append(get_interaction_point_C_d(vehicle.wake_history.history[0].rotor_wakes[rotor_idx_2].interaction_perRotor[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction[blade_idx2]))
									# l[rotor_idx][rotor_idx_2][blade_idx][blade_idx2].append(get_interaction_point_l(vehicle.wake_history.history[0].rotor_wakes[rotor_idx_2].interaction_perRotor[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction[blade_idx2]))
									# secLen[rotor_idx][rotor_idx_2][blade_idx][blade_idx2].append(get_interaction_point_secLen(vehicle.wake_history.history[0].rotor_wakes[rotor_idx_2].interaction_perRotor[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction[blade_idx2]))
									
									tvi = vehicle.wake_history.history[0].rotor_wakes[rotor_idx_2].interaction_perRotor[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction[blade_idx2]
									temp_out = []
									tmep_out2 = []
									fill_interactionPt_wake_idx(tvi, wake_idx[rotor_idx][rotor_idx_2][blade_idx][blade_idx2])

									tmep_out2 = get_interactionPt_wake_idx(tvi)

									wake_idx[rotor_idx][rotor_idx_2][blade_idx][blade_idx2].append(get_interactionPt_wake_idx(tvi))
									bladeSec_idx[rotor_idx][rotor_idx_2][blade_idx][blade_idx2].append(get_interactionPt_bladeSec_idx(vehicle.wake_history.history[0].rotor_wakes[rotor_idx_2].interaction_perRotor[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction[blade_idx2]))
									blade_directionVec[rotor_idx][rotor_idx_2][blade_idx][blade_idx2].append(get_interaction_point_r_blade(vehicle.wake_history.history[0].rotor_wakes[rotor_idx_2].interaction_perRotor[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction[blade_idx2]))
									directionVec_bv[rotor_idx][rotor_idx_2][blade_idx][blade_idx2].append(get_interaction_point_r_blade_v(vehicle.wake_history.history[0].rotor_wakes[rotor_idx_2].interaction_perRotor[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction[blade_idx2]))
									vortex_directionVec[rotor_idx][rotor_idx_2][blade_idx][blade_idx2].append(get_interaction_point_r_vortex(vehicle.wake_history.history[0].rotor_wakes[rotor_idx_2].interaction_perRotor[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction[blade_idx2]))
									normalVec[rotor_idx][rotor_idx_2][blade_idx][blade_idx2].append(get_interaction_point_blade_normal(vehicle.wake_history.history[0].rotor_wakes[rotor_idx_2].interaction_perRotor[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction[blade_idx2]))
									missDist[rotor_idx][rotor_idx_2][blade_idx][blade_idx2].append(get_interaction_point_miss_dist(vehicle.wake_history.history[0].rotor_wakes[rotor_idx_2].interaction_perRotor[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction[blade_idx2]))
									gammaSec[rotor_idx][rotor_idx_2][blade_idx][blade_idx2].append(get_interaction_point_gamma_sec(vehicle.wake_history.history[0].rotor_wakes[rotor_idx_2].interaction_perRotor[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction[blade_idx2]))
									gammaVortex[rotor_idx][rotor_idx_2][blade_idx][blade_idx2].append(get_interaction_point_gamma_w(vehicle.wake_history.history[0].rotor_wakes[rotor_idx_2].interaction_perRotor[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction[blade_idx2]))
									r_c[rotor_idx][rotor_idx_2][blade_idx][blade_idx2].append(get_interaction_point_r_c(vehicle.wake_history.history[0].rotor_wakes[rotor_idx_2].interaction_perRotor[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction[blade_idx2]))
									C_d[rotor_idx][rotor_idx_2][blade_idx][blade_idx2].append(get_interaction_point_C_d(vehicle.wake_history.history[0].rotor_wakes[rotor_idx_2].interaction_perRotor[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction[blade_idx2]))
									l[rotor_idx][rotor_idx_2][blade_idx][blade_idx2].append(get_interaction_point_l(vehicle.wake_history.history[0].rotor_wakes[rotor_idx_2].interaction_perRotor[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction[blade_idx2]))
									secLen[rotor_idx][rotor_idx_2][blade_idx][blade_idx2].append(get_interaction_point_secLen(vehicle.wake_history.history[0].rotor_wakes[rotor_idx_2].interaction_perRotor[rotor_idx].blade_vortex_interaction[blade_idx].tip_vortex_interaction[blade_idx2]))
									
														 
						fill_aoaf(blade, aoa_array)
						fill_u_tf(blade, u_t)
						fill_u_pf(blade, u_p_array)

						u = abs(omegas[rotor_idx])*radii[rotor_idx]*np.sqrt(u_t**2.0 + u_p_array**2.0)
						u = u.astype(dtype=np.single)

						append_bpm_data(bpm_files[rotor_idx][blade_idx], loading_data.time, aoa_array, 2.0*math.pi, u)

					if (trackBWIevents):

						K1 = 0.09
						K2 = 0.11

						for rotor_idx_2 in range(num_rotors):
							size2 = len(wake_idx[rotor_idx][rotor_idx_2][blade_idx][blade_idx2])
							for blade_idx2 in range(numBladesArray[rotor_idx_2]):
								for blade_idx in range(numBladesArray[rotor_idx]):
								
									if len(wake_idx[rotor_idx][rotor_idx_2][blade_idx][blade_idx2])>1:
										print('change the logic')
									size = len(wake_idx[rotor_idx][rotor_idx_2][blade_idx][blade_idx2][0])

									for i in range (0,size):
										if wake_idx[rotor_idx][rotor_idx_2][blade_idx][blade_idx2][0][i] > 2:
											interaction[blade_idx2].blade_idx.extend([blade_idx])
											interaction[blade_idx2].wake_idx.append(wake_idx[rotor_idx][rotor_idx_2][blade_idx][blade_idx2][0][i])
											interaction[blade_idx2].directionVec_bv.append(directionVec_bv[rotor_idx][rotor_idx_2][blade_idx][blade_idx2][0][i])
											interaction[blade_idx2].secIdx.append(bladeSec_idx[rotor_idx][rotor_idx_2][blade_idx][blade_idx2][0][i])
											interaction[blade_idx2].vortex_directionVec.append(vortex_directionVec[rotor_idx][rotor_idx_2][blade_idx][blade_idx2][0][i])
											interaction[blade_idx2].blade_directionVec.append(blade_directionVec[rotor_idx][rotor_idx_2][blade_idx][blade_idx2][0][i])
											interaction[blade_idx2].normalVec.append(normalVec[rotor_idx][rotor_idx_2][blade_idx][blade_idx2][0][i])
											interaction[blade_idx2].missDist.append(missDist[rotor_idx][rotor_idx_2][blade_idx][blade_idx2][0][i])
											interaction[blade_idx2].gamma.append(gammaVortex[rotor_idx][rotor_idx_2][blade_idx][blade_idx2][0][i])
											interaction[blade_idx2].Cd.append(C_d[rotor_idx][rotor_idx_2][blade_idx][blade_idx2][0][i])
											interaction[blade_idx2].r_c.append(r_c[rotor_idx][rotor_idx_2][blade_idx][blade_idx2][0][i])
											interaction[blade_idx2].l.append(l[rotor_idx][rotor_idx_2][blade_idx][blade_idx2][0][i])
										else:
											interaction[blade_idx2].gamma0 = gammaVortex[rotor_idx][rotor_idx_2][blade_idx][blade_idx2][0][i]

									
										#append_bwi_data(bwi_files[rotor_idx][blade_idx], loading_data.time, aoa_array, 2.0*math.pi, u)
									#print('2. blade_idx2:', blade_idx2, 'rotor_idx_2:', rotor_idx_2, 'blade_idx:', blade_idx, 'rotor_idx:', rotor_idx)

							
							#flat_wake_idx = [w[0] if isinstance(w, list) else w for w in interaction[blade_idx].wake_idx]
						for blade_idx in range(0,rotor.blade_states.length()):
							sorted_indices = [int(idx) for idx in np.argsort(interaction[blade_idx].wake_idx)]
							interaction[blade_idx].wake_idx = [interaction[blade_idx].wake_idx[i] for i in sorted_indices]
							interaction[blade_idx].blade_idx = [interaction[blade_idx].blade_idx[i] for i in sorted_indices]
							interaction[blade_idx].directionVec_bv = [interaction[blade_idx].directionVec_bv[i] for i in sorted_indices]
							interaction[blade_idx].secIdx = [interaction[blade_idx].secIdx[i] for i in sorted_indices]
							interaction[blade_idx].vortex_directionVec = [interaction[blade_idx].vortex_directionVec[i] for i in sorted_indices]
							interaction[blade_idx].blade_directionVec = [interaction[blade_idx].blade_directionVec[i] for i in sorted_indices]
							interaction[blade_idx].normalVec = [interaction[blade_idx].normalVec[i] for i in sorted_indices]
							interaction[blade_idx].missDist = [interaction[blade_idx].missDist[i] for i in sorted_indices]
							interaction[blade_idx].gamma = [interaction[blade_idx].gamma[i] for i in sorted_indices]
							interaction[blade_idx].Cd = [interaction[blade_idx].Cd[i] for i in sorted_indices]
							interaction[blade_idx].r_c = [interaction[blade_idx].r_c[i] for i in sorted_indices]
							interaction[blade_idx].l = [interaction[blade_idx].l[i] for i in sorted_indices]
							#interaction[blade_idx].secLen = [interaction[blade_idx].secLen[i] for i in sorted_indices]

							for i in range(0,len(interaction[blade_idx].wake_idx)):
								a = interaction[blade_idx].blade_directionVec[i]
								b = interaction[blade_idx].vortex_directionVec[i]
								c = interaction[blade_idx].normalVec[i]
								d = interaction[blade_idx].directionVec_bv[i]
								wake_idx_values = interaction[blade_idx].wake_idx[i]
								
								blade_idx_values = np.array(interaction[blade_idx].blade_idx[i])
								secIdx_values = np.array(interaction[blade_idx].secIdx[i])
								missDist_values = np.array(interaction[blade_idx].missDist[i])
								gamma_values = np.array(interaction[blade_idx].gamma[i])
								gamma0 = interaction[blade_idx].gamma0
								Cd_values = np.array(interaction[blade_idx].Cd[i])
								r_c_values = np.array(interaction[blade_idx].r_c[i])
								l_values = np.array(interaction[blade_idx].l[i])
								#secLen = np.array(interaction[tipV_idx].secLen)
								
								angle1 = np.dot(np.ravel(a), np.ravel(b)) / (np.linalg.norm(a) * np.linalg.norm(b))
								angle2 = np.dot(np.ravel(c), np.ravel(b)) / (np.linalg.norm(c) * np.linalg.norm(b))
								
								value = np.dot(np.ravel(d), np.ravel(c)) / np.linalg.norm(c)
								if abs(angle1) < 0.35 and abs(angle2) < 0.35:
									angle_offset = 360/rotor.blade_states.length()
									angle = acoustic_iteration + (180 + blade_idx_values*angle_offset)%360

									perpInteraction[blade_idx].psi.append(acoustic_iteration + angle)	
									perpInteraction[blade_idx].wake_idx.append(wake_idx_values)
									perpInteraction[blade_idx].blade_idx.append(blade_idx_values)
									perpInteraction[blade_idx].blade_directionVec.append(a)
									perpInteraction[blade_idx].secIdx.append(secIdx_values)
									perpInteraction[blade_idx].vortex_directionVec.append(b)
									perpInteraction[blade_idx].normalVec.append(c)
									perpInteraction[blade_idx].missDist.append(missDist_values)
									perpInteraction[blade_idx].missDist2.append(abs(value))
									perpInteraction[blade_idx].gamma.append(gamma_values)
									perpInteraction[blade_idx].gamma0 = gamma0
									perpInteraction[blade_idx].Cd.append(Cd_values)
									perpInteraction[blade_idx].r_c.append(r_c_values)
									perpInteraction[blade_idx].l.append(l_values)
									#perpInteraction[tipV_idx].secLength.append(np.sqrt(secLen[i]))
							
						for tipV_idx in range(0,rotor.blade_states.length()):
							#print('tipV_idx:', tipV_idx)
							if perpInteraction[tipV_idx].wake_idx:
								a = perpInteraction[tipV_idx].blade_directionVec
								b = perpInteraction[tipV_idx].vortex_directionVec
								c = perpInteraction[tipV_idx].normalVec
								wake_idx_values = np.array(perpInteraction[tipV_idx].wake_idx)
								psi_values = np.array(perpInteraction[tipV_idx].psi)
								secIdx_values = np.array(perpInteraction[tipV_idx].secIdx)
								missDist_values = np.array(perpInteraction[tipV_idx].missDist)
								missDist2 = np.array(perpInteraction[tipV_idx].missDist2)
								gamma_values = np.array(perpInteraction[tipV_idx].gamma)
								gamma0 = perpInteraction[tipV_idx].gamma0
								Cd_values = np.array(perpInteraction[tipV_idx].Cd)
								r_c_values = np.array(perpInteraction[tipV_idx].r_c)
								l_values = np.array(perpInteraction[tipV_idx].l)
								blade_idx_values = np.array(perpInteraction[tipV_idx].blade_idx)
								#secLen_values = np.array(perpInteraction[tipV_idx].secLength)
								theta = np.zeros(16)									
								for i in range(len(wake_idx_values)):
									blade_idxx = blade_idx_values[i]
									
									#perpInteraction_perBlade[blade_idxx].secLength.append(secLen_values[i])

									Vx = flight_condition["V_inf"] + omegas[rotor_idx]*radii[rotor_idx]* np.sin(np.radians(acoustic_iteration+ psi_values[i] - 1))
									Vy = omegas[rotor_idx]*radii[rotor_idx]* np.cos(np.radians(acoustic_iteration + psi_values[i] - 1))
									Uref = np.sqrt(Vx**2 + Vy**2)
									
								
									# Assuming constant chord, in the future change this to take chord length at the scetion of interaction
									theta[i] = Cd_values[i]*real_c[0]/ 2
									theta0 = 0.029*real_c[0]/ 2
									L0 = 0.32 * theta0 * np.sqrt(l_values[i] / theta0 + 380)
									perpInteraction_perBlade[blade_idxx].psi.append(psi_values[i])
									perpInteraction_perBlade[blade_idxx].ID.append(i+1)
									perpInteraction_perBlade[blade_idxx].wake_idx.append(wake_idx_values[i])
									perpInteraction_perBlade[blade_idxx].tipV_idx.append(tipV_idx)
									perpInteraction_perBlade[blade_idxx].blade_directionVec.append(np.array(a[i]))
									perpInteraction_perBlade[blade_idxx].secIdx.append(secIdx_values[i])
									perpInteraction_perBlade[blade_idxx].vortex_directionVec.append(np.array(b[i]))
									perpInteraction_perBlade[blade_idxx].normalVec.append(np.array(c[i]))
									perpInteraction_perBlade[blade_idxx].missDist.append(missDist_values[i])
									perpInteraction_perBlade[blade_idxx].missDist2.append(missDist2[i])
									perpInteraction_perBlade[blade_idxx].gamma.append(gamma_values[i])
									perpInteraction_perBlade[blade_idxx].gamma00.append(gamma0)
									perpInteraction_perBlade[blade_idxx].Cd.append(Cd_values[i])
									perpInteraction_perBlade[blade_idxx].r_c.append(r_c_values[i])
									perpInteraction_perBlade[blade_idxx].l.append(l_values[i])
									perpInteraction_perBlade[blade_idxx].Uref.append(Uref)
									perpInteraction_perBlade[blade_idxx].L0.append(L0)
										
									if i > 0:
										L0_1 = 0.32 * theta[0] * np.sqrt(l_values[0] / theta[0] + 380)
										L0_2 = 0.32 * theta[1] * np.sqrt(l_values[1] / theta[1] + 380)
										value = K1 * np.sqrt(l_values[0]/real_c[0]) + abs(K2*gamma0)/(Uref*real_c[0])*((l_values[0] * L0_1)/real_c[0]**2 + (l_values[1]*L0_2)/real_c[0]**2)
										perpInteraction_perBlade[blade_idxx].a0.append(value)
										if (value**2 - missDist2[i]**2) > 0.0:
											perpInteraction_perBlade[blade_idxx].b_e.append(np.sqrt(value**2 - missDist2[i]**2))
										else:
											perpInteraction_perBlade[blade_idxx].b_e.append(0.0)
									else:
										perpInteraction_perBlade[blade_idxx].a0.append(r_c_values[i]) 
										perpInteraction_perBlade[blade_idxx].b_e.append(r_c_values[i])
										
										# value = r_c_values[i]**2 - missDist2[i]**2
										# if value > 0:
										# 	perpInteraction_perBlade[blade_idxx].b_e.append(value)
										# else:
										# 	value = 0.0
										# 	perpInteraction_perBlade[blade_idxx].b_e.append(r_c_values[i])
									
						
						for blade_idxx in range(0,rotor.blade_states.length()): 
							
							ID = []
							psi = []
							sec_idx = []
							wake_angle = []
							missDist2 = []
							vortex_len = []
							L0 = []
							b_e = []
							Uref = []
							
							ID_full = perpInteraction_perBlade[blade_idxx].ID
							psi_full = perpInteraction_perBlade[blade_idxx].psi
							sec_idx_full = perpInteraction_perBlade[blade_idxx].secIdx
							wake_angle_full = perpInteraction_perBlade[blade_idxx].wake_idx
							missDist2_full = perpInteraction_perBlade[blade_idxx].missDist2
							vortex_len_full = perpInteraction_perBlade[blade_idxx].l
							L0_full = perpInteraction_perBlade[blade_idxx].L0
							b_e_full = perpInteraction_perBlade[blade_idxx].b_e
							Uref_full = perpInteraction_perBlade[blade_idxx].Uref
							
							filtered_indices = [i for i, id_val in enumerate(ID_full) if id_val != 1]
							nInteractions = len(filtered_indices)

							if nInteractions > 0:
								ID = [ID_full[i] for i in filtered_indices]
								psi = [psi_full[i] for i in filtered_indices]
								sec_idx = [sec_idx_full[i] for i in filtered_indices]
								wake_angle = [wake_angle_full[i] for i in filtered_indices]
								missDist2 = [missDist2_full[i] for i in filtered_indices]
								vortex_len = [vortex_len_full[i] for i in filtered_indices]
								L0 = [L0_full[i] for i in filtered_indices]
								b_e = [b_e_full[i] for i in filtered_indices]
								Uref = [Uref_full[i] for i in filtered_indices]
								#print('b_e:', b_e)
							# print('writing rotor_idx:', rotor_idx, 'blade_idxx:', blade_idxx)
							append_bwi_data(bwi_files[rotor_idx][blade_idxx], loading_data.time, nInteractions, ID, psi, sec_idx, wake_angle, missDist2, vortex_len, L0, b_e, Uref)		
							

				spanwise_element_iteration = spanwise_element_iteration + 1
				
				if track_wake_element:
					for t_idx in range(len(target_y_slices)):
						for rotor_idx in range(num_rotors):
							if wake_element_found[rotor_idx][t_idx]:
								fill_wake_x_component(vehicle.wake_history.history[0].rotor_wakes[0].tip_vortices[wake_element_blade[rotor_idx, t_idx]], temp_wake_array_x)
								fill_wake_z_component(vehicle.wake_history.history[0].rotor_wakes[0].tip_vortices[wake_element_blade[rotor_idx, t_idx]], temp_wake_array_z)
								fill_wake_r_c_component(vehicle.wake_history.history[0].rotor_wakes[0].tip_vortices[wake_element_blade[rotor_idx, t_idx]], temp_wake_array_r_c)

								if wake_element_index[rotor_idx, t_idx] < min(int(post_conv_revolutions*iter_per_rev), wake_lengths[0]):
									wake_element_trajectory[rotor_idx, t_idx, 0, wake_element_index[rotor_idx, t_idx]] = temp_wake_array_x[wake_element_index[rotor_idx, t_idx]]/radii[rotor_idx]
									wake_element_trajectory[rotor_idx, t_idx, 1, wake_element_index[rotor_idx, t_idx]] = temp_wake_array_z[wake_element_index[rotor_idx, t_idx]]/radii[rotor_idx]
									wake_element_core_size[rotor_idx, t_idx, wake_element_index[rotor_idx, t_idx]] = temp_wake_array_r_c[wake_element_index[rotor_idx, t_idx]]

									wake_element_index[rotor_idx, t_idx] = wake_element_index[rotor_idx, t_idx] + 1

				if track_piv_window:
					for rotor_idx in range(num_rotors):
						for blade_idx in range(num_blades[rotor_idx]):
							fill_wake_x_component(vehicle.wake_history.history[0].rotor_wakes[rotor_idx].tip_vortices[blade_idx], temp_wake_array_x)
							fill_wake_y_component(vehicle.wake_history.history[0].rotor_wakes[rotor_idx].tip_vortices[blade_idx], temp_wake_array_y)
							fill_wake_z_component(vehicle.wake_history.history[0].rotor_wakes[rotor_idx].tip_vortices[blade_idx], temp_wake_array_z)

							for t_idx in range(len(piv_slices)):

								in_y_plane = []
								in_x_left = []
								in_x_right = []
								in_x_bounds = []
								in_z_bottom = []
								in_z_top = []
								in_z_bounds = []
								in_piv_window = []

								in_y_plane = np.abs(temp_wake_array_y - piv_slices[t_idx]) < 0.05

								in_x_left = temp_wake_array_x > piv_window_x[0]
								in_x_right = temp_wake_array_x < piv_window_x[1]
								in_x_bounds = np.logical_and(in_x_left, in_x_right)

								in_z_bottom = temp_wake_array_z > piv_window_z[0]
								in_z_top = temp_wake_array_z < piv_window_z[1]
								in_z_bounds = np.logical_and(in_z_bottom, in_z_top)

								in_piv_window = np.logical_and(in_x_bounds, in_y_plane, in_z_bounds)

								if np.any(in_piv_window):
									global_pos = Vec3([temp_wake_array_x[in_piv_window][0], temp_wake_array_y[in_piv_window][0], temp_wake_array_z[in_piv_window][0]])
									rotor_local_pos = vehicle.aircraft.rotors[results["piv_window"]["rotor"]].frame.parent.global_to_local(global_pos)

									wake_element_piv[rotor_idx, blade_idx, t_idx, 0, piv_window_index[rotor_idx, blade_idx, t_idx]] = rotor_local_pos[0]
									wake_element_piv[rotor_idx, blade_idx, t_idx, 1, piv_window_index[rotor_idx, blade_idx, t_idx]] = rotor_local_pos[1]

									piv_window_index[rotor_idx, blade_idx, t_idx] = piv_window_index[rotor_idx, blade_idx, t_idx] + 1

				acoustic_iteration = acoustic_iteration + 1

			iteration = iteration + 1
			print('iteration:', iteration, 'acoustic_iteration:', acoustic_iteration)
		for rotor_idx in range(num_rotors):
			for blade_idx in range(num_blades[rotor_idx]):
				close_loading_file(loading_files[rotor_idx][blade_idx])
				if(trackBWIevents):
					close_bwi_file(bwi_files[rotor_idx][blade_idx])

	log_file.write("Sim done\n")

	result_dictionary = {}

	if not do_compute:
		try:
			result_dictionary = scipy.io.loadmat(f"{output_base}/results.mat")
			omegas = result_dictionary["omegas"][0]
			dt = result_dictionary["dt"][0][0]
			rotor_phases = result_dictionary["rotor_phases"][0]
			print(f'{rotor_phases}')
			print(f'{omegas}')
			print(f'{dt}')
			for rotor_idx in range(vehicle.aircraft.rotors.length()):
				vehicle.input_state.rotor_inputs[rotor_idx].angular_velocity = omegas[rotor_idx]

		except:
			raise FileExistsError("Failed to find results.mat from previous run")

	# Nitya, 09.16
	# if trackBWIevents:
	# 	result_dictionary['wake_idx'] = np.asarray(wake_idx, dtype=object)
	# 	result_dictionary['bladeSec_idx'] = np.asarray(bladeSec_idx, dtype=object)
	# 	result_dictionary['wake_miss_dist'] = np.asarray(missDist, dtype=object)
	# 	result_dictionary['r_c'] = np.asarray(r_c, dtype=object)
	# 	result_dictionary['gamma_sec'] = np.asarray(gammaSec, dtype=object)
	# 	result_dictionary['gamma_vortex'] = np.asarray(gammaVortex, dtype=object)
	# 	result_dictionary['C_d'] = np.asarray(C_d, dtype=object)
	# 	result_dictionary['l'] = np.asarray(l, dtype=object)
	# 	result_dictionary['secLen'] = np.asarray(secLen, dtype=object)
	# 	result_dictionary['blade_directionVec'] = np.asarray(blade_directionVec, dtype=object)
	# 	result_dictionary['vortex_directionVec'] = np.asarray(vortex_directionVec, dtype=object)
	# 	result_dictionary['directionVec_bv'] = np.asarray(directionVec_bv, dtype=object)
	# 	result_dictionary['normalVec'] = np.asarray(normalVec, dtype=object)

	if elastic_twist is not None:
		result_dictionary['elastic_twist_array'] = elastic_twist_array
		result_dictionary['blade_twist_array'] = blade_twist_array
		result_dictionary['blade_twist_azimuth'] = blade_twist_azimuth

		result_dictionary['collective_pitch_array'] = collective_pitch_array
		result_dictionary['sin_pitch_array'] = sin_pitch_array
		result_dictionary['cos_pitch_array'] = cos_pitch_array
		result_dictionary['hhc_pitch_array'] = hhc_pitch_array

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
			result_dictionary['span_element_af_loading'] = span_element_af_loading
			result_dictionary['span_element_loading'] = span_element_loading
			result_dictionary['span_element_aoa_eff'] = span_element_aoa_eff
			result_dictionary['span_element_aoa'] = span_element_aoa
			result_dictionary['span_element_up'] = span_element_up
			result_dictionary['span_element_inflow_angle'] = span_element_inflow_angle
			result_dictionary['span_element_theta'] = span_element_theta
			result_dictionary['span_element_gamma'] = span_element_gamma
			result_dictionary['span_element_azimuth'] = span_element_azimuth

		if track_piv_window:
			result_dictionary['wake_element_piv'] = wake_element_piv

		result_dictionary["omegas"] = omegas
		result_dictionary["dt"] = dt
		result_dictionary["average_powers"] = average_Qs
		result_dictionary["average_torques"] = [average_Qs[rotor_idx]/abs(omegas[rotor_idx]) for rotor_idx in range(num_rotors)]
		result_dictionary["rotor_phases"] = rotor_phases

		result_dictionary["blade_inflow_distribution"] = blade_inflow_distribution
		result_dictionary["blade_loading_distribution"] = blade_loading_distribution

		result_dictionary["blade_induced_drag_distribution"] = blade_induced_drag_distribution
		result_dictionary["blade_profile_drag_distribution"] = blade_profile_drag_distribution
		result_dictionary["blade_dynamic_aoa_distribution"] = blade_dynamic_aoa_distribution

	namelists = []

	if (acoustics is not None) and (observer is not None):

		min_omega = np.min(np.abs(omegas))

		tau_min = 0.5*(1.0*math.pi/min_omega)
		tau_max = tau_min + (post_conv_revolutions - 0.5)*(2.0*math.pi/min_omega)

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

		elif observer["type"] == "external_file":
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

		for rotor_idx, rotor in enumerate(vehicle.aircraft.rotors):
			wopwop_motion[rotor.frame.name]["omega"] = vehicle.input_state.rotor_inputs[rotor_idx].angular_velocity

		for rotor_idx in range(num_rotors):
			wopwop_case_path = f'{wopwop_output_path}/rotor_{rotor_idx}/'

			if not path.isdir(wopwop_case_path):
				makedirs(wopwop_case_path, exist_ok=True)

			#def generate_wopwop_namelist(atmo, dt, V_inf, iterations, aoa, t_min, t_max, nt, observer_config, acoustics_config, wopwop_data_path, sos, aircraft, rotors, wopwop_motion, ac_input, wopwop_case_path, rotor_phases):			
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
				wopwop_case_path,
				[rotor_phases[rotor_idx]],
				geom_directory
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
			wopwop_case_path,
			rotor_phases,
			geom_directory
		)

		namelists.append(namelist)

	return average_C_Ts, namelists, result_dictionary
