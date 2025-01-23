#!/usr/bin/env python3

import sys
import os

sys.path.insert(0, f'{os.path.dirname(os.path.realpath(__file__))}/../../dependencies/wopwopd')

from libwopwopd import *
import json
import numpy as np
import math

import matplotlib.pyplot as plt

import scipy.io as sio
import scipy.fft as sft
import scipy.signal as sig

import pyjson5

import matplotlib.font_manager as font_manager
import matplotlib

font_dir = '/usr/share/fonts/truetype/msttcorefonts/times.ttf'
font_manager.fontManager.addfont(font_dir)
prop = font_manager.FontProperties(fname=font_dir)

font_size0 = 18
font_size1 = 17
font_size2 = 16
font_size3 = 12
font_size4 = 5

label_font = {'fontname':'Times New Roman', 'size':f'{font_size1}', 'color':'black', 'weight':'normal',
              'verticalalignment':'bottom'}

title_font = {'fontname':'Times New Roman', 'size':f'{font_size1}', 'color':'black', 'weight':'normal',
              'verticalalignment':'bottom'}

title_font2 = {'fontname':'Times New Roman', 'size':f'{font_size2}', 'color':'black', 'weight':'normal',
              'verticalalignment':'bottom'}

legend_font = {'fontname':'Times New Roman', 'size':f'{font_size4}', 'color':'black', 'weight':'normal',
              'verticalalignment':'bottom'}

RADPS_2_HZ = 0.1591549433
R = 2

TITLE_DICT = {
	"BL": "Baseline",
	"MN": "Minimum Noise",
	"MV": "Minimum Vibration"
}

def loadmat(filename: str):
	'''
	this function should be called instead of direct spio.loadmat
	as it cures the problem of not properly recovering python dictionaries
	from mat files. It calls the function check keys to cure all entries
	which are still mat-objects
	'''
	data = sio.loadmat(filename, struct_as_record=False, squeeze_me=True)
	return _check_keys(data)

def _check_keys(dict: dict):
	'''
	checks if entries in dictionary are mat-objects. If yes
	todict is called to change them to nested dictionaries
	'''
	for key in dict:
		if isinstance(dict[key], sio.matlab.mat_struct):
			dict[key] = _todict(dict[key])
	return dict        

def _todict(matobj):
	'''
	A recursive function which constructs from matobjects nested dictionaries
	'''
	dict = {}
	for strg in matobj._fieldnames:
		elem = matobj.__dict__[strg]
		if isinstance(elem, sio.matlab.mat_struct):
			dict[strg] = _todict(elem)
		else:
			dict[strg] = elem
	return dict

def read_hart_contour_tecplot(filename: str):
	with open(filename, 'r') as tecfile:
		_ = tecfile.readline()
		_ = tecfile.readline()
		
		zone_header = tecfile.readline()
		zone_dictionary = {s[0].strip(): ''.join(s[1:]) for s in [s.split("=") for s in [s.strip() for s in zone_header[4:].split(',')]]}

		i_max = int(zone_dictionary["I"])
		j_max = int(zone_dictionary["J"])

		x_grid = np.zeros((i_max, j_max))
		y_grid = np.zeros((i_max, j_max))
		mid_freq = np.zeros((i_max, j_max))

		#i = 0
		#j = 0
		for j in range(j_max):
			for i in range(i_max):
				str_vals = tecfile.readline().split()
				x_grid[i, j] = float(str_vals[0])/R
				y_grid[i, j] = -float(str_vals[1])/R
				mid_freq[i, j] = float(str_vals[3])
		
		return x_grid, y_grid, mid_freq

def read_hart_wake_tecplot(filename: str):
	print("reading hart wake")

	with open(filename, 'r') as tecfile:
		_ = tecfile.readline()
		_ = tecfile.readline()

		wake_points = []

		done = False
		while not done:
			zone_header = tecfile.readline().strip()
			if "ZONE" in zone_header[0:4]:
				zone_dictionary = {s[0]: ''.join(s[1:]) for s in [s.split("=") for s in [s.strip() for s in zone_header[4:].split(',')]]}
				num_points = int(zone_dictionary["I"])
				# Skip root vorticies
				if "inboard" not in zone_dictionary["T"]:
					for _ in range(num_points):
						point_line = tecfile.readline().split()
						data = (
							int(float(point_line[1])),
							float(point_line[3]),
							float(point_line[5])
						)
						wake_points.append(data)
				else:
					for _ in range(num_points):
						_ = tecfile.readline().split()
			else:
				done = True

	return wake_points

def read_hart_blade_data_tecplot(filename: str):
	with open(filename, 'r') as tecfile:
		_ = tecfile.readline()
		_ = tecfile.readline()

		zone_header = tecfile.readline().strip()
		zone_dictionary = {s[0]: ''.join(s[1:]) for s in [s.split("=") for s in [s.strip() for s in zone_header[4:].split(',')]]}

		num_points = int(zone_dictionary["I"])
		blade_data = np.zeros(num_points)

		for p_idx in range(num_points):
			blade_data[p_idx] = float(tecfile.readline().strip())

	return blade_data

cmap_lines = matplotlib.colors.LinearSegmentedColormap.from_list("", ["black", "black"])

def plot_acoustic_contours(plot_name: str):
	x_grid, y_grid, measured = read_hart_contour_tecplot(f'{os.path.dirname(os.path.realpath(__file__))}/{plot_name.lower()}-contour-meas.tec')

	wopwop_results = parse_wopwop_results(f'{os.path.dirname(os.path.realpath(__file__))}/{plot_name.upper()}/acoustics/full_system', 'case.nam')

	i_max = len(wopwop_results.oaspl_db_grid.obs_x)
	j_max = len(wopwop_results.oaspl_db_grid.obs_x[0])

	print(f'i_max: {i_max}')
	print(f'j_max: {j_max}')

	oaspl_linear = [oaspl_db.functions[2].data[0] for oaspl_db in wopwop_results.oaspl_db]

	oaspl_db = [[oaspl_linear[i*j_max + j] for i in range(i_max)] for j in range(j_max)]

	#print(wopwop_results.oaspl_db_grid.obs_y[0])
	#print(wopwop_results.oaspl_db_grid.obs_x[1])

	x = wopwop_results.oaspl_db_grid.obs_x[0]
	y = [_y[0] for _y in wopwop_results.oaspl_db_grid.obs_y]
	y.reverse()

	#print(x)
	offset = x[0] - -4.0

	clevels = np.linspace(85, 119, 18)

	#light_rainbow = cmap_map(lambda x: x/2 + 0.5, matplotlib.cm.rainbow)

	print([(_x - offset)/R for _x in x])
	print([_y/R for _y in y])
	print(x_grid[0,:])
	print(y_grid[:,0])
	#print(x)

	fig = plt.figure()
	ax0 = plt.subplot(121)
	plt2 = plt.contour([-_y/R for _y in y], [-(_x - offset)/R for _x in x], oaspl_db, levels=clevels, cmap=cmap_lines, linewidths=0.5)
	plt1 = plt.contourf([-_y/R for _y in y], [-(_x - offset)/R for _x in x], oaspl_db, levels=clevels, cmap='Blues')
	#plt.clabel(plt2, clevels, inline=True, colors='k', fontsize=5)
	clabels = plt.clabel(plt2, clevels, colors='k', fontsize=font_size4)
	plt.ylabel('$x/R$', **label_font)
	plt.xlabel('$y/R$', labelpad=20, **label_font)
	plt.title('Prediction', **title_font)
	#plt.axis('scaled')
	
	#plt.axis('equal')
	
	for label in clabels:
		label.set_bbox(dict(facecolor='white', edgecolor='white', pad=0.01))
	
	ax = plt.gca()
	for label in (ax.get_xticklabels() + ax.get_yticklabels()):
		label.set_fontname('Times New Roman')
		label.set_fontsize(font_size3)

	ax1 = plt.subplot(122)
	ax1.set_yticklabels([])
	#ax.set_xticklabels([])
	plt2 = plt.contour(x_grid[0,:], y_grid[:,0], measured, levels=clevels, cmap=cmap_lines, linewidths=0.5)
	plt.contourf(x_grid[0,:], y_grid[:,0], measured, levels=clevels, cmap='Blues')
	clabels = plt.clabel(plt2, clevels, colors='k', fontsize=font_size4)
	plt.title('Measured', **title_font)
	#plt.clabel(CS, clevels, inline=True)
	#plt.contourf(y_grid[:,0], x_grid[0,:], measured, levels=clevels)
	#plt.axis('equal')
	#plt.colorbar()
	plt.xlabel('$y/R$', labelpad=20, **label_font)
	#plt.ylabel('-x/R')

	for label in clabels:
		label.set_bbox(dict(facecolor='white', edgecolor='white', pad=0.01))

	ax = plt.gca()
	for label in (ax.get_xticklabels() + ax.get_yticklabels()):
		label.set_fontname('Times New Roman')
		label.set_fontsize(font_size3)

	#ax = plt.subplot(133)
	cax = fig.add_axes([0.93, 0.11, 0.02, 0.77])
	cb = fig.colorbar(plt1, cax, orientation='vertical')#, label='BVI SPL [dB]', **label_font)
	cb.set_label('BVI SPL [dB]', labelpad=20, **label_font)

	#ax = plt.gca()
	for label in (cax.get_xticklabels() + cax.get_yticklabels()):
			label.set_fontname('Times New Roman')
			label.set_fontsize(font_size3)

	plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/Hart_{plot_name.upper()}.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
	#plt.show()

def plot_blade_normal_pressures(plot_name: str):

	with open(f'{os.path.dirname(os.path.realpath(__file__))}/../hart_ii_params.json5') as param_file:
		params = pyjson5.load(param_file)

	flight_condition = list(filter(lambda x: x["name"] == plot_name, params['flight_conditions']))[0]

	#omega = flight_condition["omegas"][0]
	omegas = []
	for m in flight_condition["motion"]:
		if ('blade_element_func' in m) and (m['blade_element_func'] == "rotation"):
			omegas.append(m['omega'])

	omega = omegas[0]
	
	blade_results = loadmat(f'{os.path.dirname(os.path.realpath(__file__))}/{plot_name.upper()}/results.mat')

	measured_blade_loading = read_hart_blade_data_tecplot(f'{os.path.dirname(os.path.realpath(__file__))}/CnM2_{plot_name}_ca.tec')

	measured_mean = measured_blade_loading.mean()
	measured_blade_loading = measured_blade_loading - measured_mean
	
	blade_loading = blade_results['span_element_loading']

	computed_mean = blade_loading.mean()

	blade_loading = blade_loading - computed_mean

	measured_azimuth = np.linspace(0, 360, measured_blade_loading.size)
	computational_azimuth = blade_results['span_element_azimuth']
	#computational_azimuth = np.linspace(0, 360, blade_loading.size)

	plt.figure(num=1)
	plt.plot(measured_azimuth, measured_blade_loading, 'r.-', computational_azimuth, blade_loading, 'b', linewidth=0.5, markersize=0.7)
	plt.title(f'{TITLE_DICT[plot_name.upper()]} 87% Span blade loading. Mean removed.')
	plt.legend(['Measured', 'OpenCOPTER'])
	plt.xlabel("$\psi$")
	plt.ylabel("$C_nM^2$")
	plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/Hart_{plot_name}_normal_loading.png', dpi=500, bbox_inches="tight", pad_inches=0.01)
	plt.cla()
	plt.clf()

	# 10/rev high pass
	measured_dt = (2.0*math.pi/omega)/measured_blade_loading.shape[0]
	computed_dt = (2.0*math.pi/omega)/blade_loading.shape[0]

	print(f'measured_dt: {measured_dt}')
	print(f'computed_dt: {computed_dt}')

	fblade_loading = sft.fft(blade_loading)
	fblade_freqz = sft.fftfreq(fblade_loading.shape[0], d=computed_dt)

	fmeasured_blade_loading = sft.fft(measured_blade_loading)
	fmeasured_blade_freqz = sft.fftfreq(fmeasured_blade_loading.shape[0], d=measured_dt)

	ten_per_rev = 10*omega*RADPS_2_HZ

	print(f'fblade_freqz: {fblade_freqz}')
	print(f'fmeasured_blade_freqz: {fmeasured_blade_freqz}')

	measured_sos = sig.butter(13, ten_per_rev, 'highpass', output='sos', fs=1/measured_dt)
	computed_sos = sig.butter(13, ten_per_rev, 'highpass', output='sos', fs=1/computed_dt)

	blade_loading = sig.sosfilt(computed_sos, blade_loading)
	measured_blade_loading = sig.sosfilt(measured_sos, measured_blade_loading)

	plt.figure(num=1)
	plt.plot(measured_azimuth, measured_blade_loading, 'r.-', computational_azimuth, blade_loading, 'b', linewidth=0.5, markersize=0.7)
	plt.xlabel("$\psi$")
	plt.ylabel("$C_nM^2$")
	plt.title(f'{TITLE_DICT[plot_name.upper()]}: 87% Span blade loading. 10/rev highpassed.')
	plt.legend(['Measured', 'OpenCOPTER'])
	#plt.xlim([20, 90])
	plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/Hart_{plot_name}_normal_loading_hp.png', dpi=500, bbox_inches="tight", pad_inches=0.01)
	plt.cla()
	plt.clf()

	
	aoa_eff = blade_results['span_element_aoa_eff']
	
	plt.figure(num=1)
	plt.plot(computational_azimuth, aoa_eff, linewidth=0.5)
	plt.title(f'{TITLE_DICT[plot_name.upper()]} 87% Effective AOA.')
	#plt.legend(['Measured', 'OpenCOPTER'])
	plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/Hart_{plot_name}_aoa_eff.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
	plt.cla()
	plt.clf()

	u_p = blade_results['span_element_up']
	
	plt.figure(num=1)
	plt.plot(computational_azimuth, u_p, linewidth=0.5)
	plt.title(f'{TITLE_DICT[plot_name.upper()]} 87% U_p.')
	#plt.legend(['Measured', 'OpenCOPTER'])
	plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/Hart_{plot_name}_u_p.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
	plt.cla()
	plt.clf()

	aoa = blade_results['span_element_aoa']
	
	plt.figure(num=1)
	plt.plot(computational_azimuth, aoa, linewidth=0.5)
	plt.title(f'{TITLE_DICT[plot_name.upper()]} 87% AOA.')
	#plt.legend(['Measured', 'OpenCOPTER'])
	plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/Hart_{plot_name}_aoa.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
	plt.cla()
	plt.clf()

	theta = blade_results['span_element_theta']
	
	plt.figure(num=1)
	plt.plot(computational_azimuth, theta, linewidth=0.5)
	plt.title(f'{TITLE_DICT[plot_name.upper()]} 87% theta.')
	#plt.legend(['Measured', 'OpenCOPTER'])
	plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/Hart_{plot_name}_theta.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
	plt.cla()
	plt.clf()

	inflow_angle = blade_results['span_element_inflow_angle']
	
	plt.figure(num=1)
	plt.plot(computational_azimuth, inflow_angle, linewidth=0.5)
	plt.title(f'{TITLE_DICT[plot_name.upper()]} 87% inflow angle.')
	#plt.legend(['Measured', 'OpenCOPTER'])
	plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/Hart_{plot_name}_inflow_angle.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
	plt.cla()
	plt.clf()

	gamma = blade_results['span_element_gamma']
	
	plt.figure(num=1)
	plt.plot(computational_azimuth, gamma, linewidth=0.5)
	plt.title(f'{TITLE_DICT[plot_name.upper()]} 87% gamma.')
	#plt.legend(['Measured', 'OpenCOPTER'])
	plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/Hart_{plot_name}_gamma.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
	plt.cla()
	plt.clf()

def plot_wake_trajectory(plot_name: str):

	# slice_indides = {
	# 	-0.4: (0, 8),
	# 	-0.55: (8, 16),
	# 	-0.7: (16, 23),
	# 	-0.85: (23, 28),
	# 	-0.97: (28, 30),
	# 	0.4: (30, 36),
	# 	0.55: (36, 42),
	# 	0.7: (42, 47),
	# 	0.85: (47, 51),
	# 	0.97: (51, 53),
	# }
	
	slice_indides = {
		0.4: (0, 8),
		0.55: (8, 16),
		0.7: (16, 23),
		0.85: (23, 28),
		0.97: (28, 30),
		-0.4: (30, 36),
		-0.55: (36, 42),
		-0.7: (42, 47),
		-0.85: (47, 51),
		-0.97: (51, 53),
	}
	print(f'Plotting wake trajectory for {plot_name}')
	wake_results = loadmat(f'{os.path.dirname(os.path.realpath(__file__))}/{plot_name.upper()}/results.mat')

	measured_wake_ret = read_hart_wake_tecplot(f'{os.path.dirname(os.path.realpath(__file__))}/dnw_ret_{plot_name.lower()}53.tec')
	measured_wake_adv = read_hart_wake_tecplot(f'{os.path.dirname(os.path.realpath(__file__))}/dnw_adv_{plot_name.lower()}53.tec')

	measured_wake = measured_wake_adv + measured_wake_ret

	measured_wake.sort(key=lambda x: x[0])

	print(f'measured_wake: {measured_wake}')
	target_y_slices = wake_results['target_y_slices']

	wake_element_index = wake_results['wake_element_index']
	wake_element_trajectory = wake_results['wake_element_trajectory']
	wake_element_core_size = wake_results['wake_element_core_size']
	
	for idx, target_y_slice in enumerate(target_y_slices):

		measured_indices = slice_indides[target_y_slice]

		wake_data_points = list(filter(lambda x: x[0] >= (measured_indices[0] + 1) and x[0] <= measured_indices[1], measured_wake))

		measured_x = np.asarray([p[1]/R for p in wake_data_points])
		measured_z = np.asarray([p[2]/R for p in wake_data_points])

		
		print(f"measured_x: {measured_x}")
		print(f"measured_z: {measured_z}")
		wake_x = wake_element_trajectory[idx, 0, 0:wake_element_index[idx]]
		wake_z = wake_element_trajectory[idx, 1, 0:wake_element_index[idx]]

		i_cos_aoa = math.cos(5.3*(math.pi/180.0))
		i_sin_aoa = math.sin(5.3*(math.pi/180.0))

		x_tppc = (wake_x*i_cos_aoa - wake_z*i_sin_aoa)
		z_tppc = -(-wake_x*i_sin_aoa - wake_z*i_cos_aoa)

		if z_tppc.size > 0 and measured_z.size > 0:
			min_idx = np.argmin(np.abs(measured_x[0] - x_tppc))
			delta = measured_z[0] - z_tppc[min_idx]

			measured_z = measured_z - delta

		plt.figure(num=1)
		plt.plot(measured_x, measured_z, 'r*', x_tppc, z_tppc, 'b-', linewidth=0.5)
		plt.title(f'{TITLE_DICT[plot_name.upper()]}: Y slice: {target_y_slice:.2f}')
		plt.legend(['Measured', 'OpenCOPTER'])
		#plt.axis('scaled')
		plt.xlim([-1.1, 1.1])
		plt.ylim([-0.1, 0.15])
		plt.xlabel("$x/R$")
		plt.ylabel("$z/R$")
		#plt.gca().invert_xaxis()
		plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/Hart_{plot_name}_wake_y{target_y_slice:.2f}.png', dpi=500, bbox_inches="tight", pad_inches=0.01)
		plt.cla()
		plt.clf()

		core_size = wake_element_core_size[idx, 0:wake_element_index[idx]]
		plt.figure(num=1)
		plt.plot(core_size, linewidth=0.5)
		plt.title(f'{TITLE_DICT[plot_name.upper()]}: Core size, Y slice: {target_y_slice:.2f}')
		#plt.legend(['Measured', 'OpenCOPTER'])
		#plt.axis('scaled')
		#plt.xlim([-1.1, 1.1])
		#plt.ylim([-0.1, 0.15])
		#plt.gca().invert_xaxis()
		plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/Hart_{plot_name}_core_size_y{target_y_slice:.2f}.png', dpi=500, bbox_inches="tight", pad_inches=0.01)
		plt.cla()
		plt.clf()

	wake_fil_core_size = wake_results['wake_0_core_size'][2,:]
	#core_size = wake_element_core_size[idx, 0:wake_element_index[idx]]
	plt.figure(num=1)
	plt.plot(wake_fil_core_size, linewidth=0.5)
	plt.title(f'{TITLE_DICT[plot_name.upper()]}: Core size, Y slice: {target_y_slice:.2f}')
	#plt.legend(['Measured', 'OpenCOPTER'])
	#plt.axis('scaled')
	#plt.xlim([-1.1, 1.1])
	#plt.ylim([-0.1, 0.15])
	#plt.gca().invert_xaxis()
	plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/Hart_{plot_name}_fil_core_size.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
	plt.cla()
	plt.clf()

def plot_blade_twist(plot_name: str):

	with open(f'{os.path.dirname(os.path.realpath(__file__))}/../hart_ii_params.json5') as param_file:
		params = pyjson5.load(param_file)

	flight_condition = list(filter(lambda x: x["name"] == plot_name, params['flight_conditions']))[0]

	omegas = []
	for m in flight_condition["motion"]:
		if ('blade_element_func' in m) and (m['blade_element_func'] == "rotation"):
			omegas.append(m['omega'])

	omega = omegas[0]
	ten_per_rev = 10*omega*RADPS_2_HZ

	blade_results = loadmat(f'{os.path.dirname(os.path.realpath(__file__))}/{plot_name.upper()}/results.mat')

	blade_twist_array = blade_results['blade_twist_array'][2,:]
	blade_twist_azimuth = blade_results['blade_twist_azimuth'][2,:]
	

	computed_dt = (2.0*math.pi/omega)/blade_twist_array.shape[0]
	computed_sos = sig.butter(13, ten_per_rev, 'highpass', output='sos', fs=1/computed_dt)

	#print(f'blade_twist_array: {blade_twist_array}')
	plt.figure(num=1)
	plt.plot(blade_twist_azimuth, blade_twist_array, linewidth=0.5)
	plt.title(f'{TITLE_DICT[plot_name.upper()]}: Blade twist over time.')
	plt.xlim([0, 360])
	plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/Hart_{plot_name}_blade_twist.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
	plt.cla()
	plt.clf()
	
	# elastic_twist_array = blade_results['elastic_twist_array'][2,:]
	
	# #print(f'blade_twist_array: {blade_twist_array}')
	# plt.figure(num=1)
	# plt.plot(blade_twist_azimuth, elastic_twist_array, linewidth=0.5)
	# plt.title(f'{TITLE_DICT[plot_name.upper()]}: Elastic twist over time.')
	# plt.xlim([0, 360])
	# plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/Hart_{plot_name}_elastic_twist.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
	# plt.cla()
	# plt.clf()

	collective_pitch_array = blade_results['collective_pitch_array'][:]
	
	#print(f'blade_twist_array: {blade_twist_array}')
	plt.figure(num=1)
	plt.plot(blade_twist_azimuth, collective_pitch_array, linewidth=0.5)
	plt.title(f'{TITLE_DICT[plot_name.upper()]}: Collective pitch over time.')
	plt.xlim([0, 360])
	plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/Hart_{plot_name}_collective.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
	plt.cla()
	plt.clf()

	sin_pitch_array = blade_results['sin_pitch_array'][2,:]
	sin_pitch_array_hp = sig.sosfilt(computed_sos, sin_pitch_array)

	#print(f'blade_twist_array: {blade_twist_array}')
	plt.figure(num=1)
	plt.plot(blade_twist_azimuth, sin_pitch_array, linewidth=0.5)
	plt.title(f'{TITLE_DICT[plot_name.upper()]}: 1S over time.')
	plt.xlim([0, 360])
	plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/Hart_{plot_name}_1s.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
	plt.cla()
	plt.clf()

	#print(f'blade_twist_array: {blade_twist_array}')
	plt.figure(num=1)
	plt.plot(blade_twist_azimuth, sin_pitch_array_hp, linewidth=0.5)
	plt.title(f'{TITLE_DICT[plot_name.upper()]}: 1S over time HP.')
	plt.xlim([0, 360])
	plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/Hart_{plot_name}_1s_hp.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
	plt.cla()
	plt.clf()

	cos_pitch_array = blade_results['cos_pitch_array'][2,:]
	cos_pitch_array_hp = sig.sosfilt(computed_sos, cos_pitch_array)

	#print(f'blade_twist_array: {blade_twist_array}')
	plt.figure(num=1)
	plt.plot(blade_twist_azimuth, cos_pitch_array, linewidth=0.5)
	plt.title(f'{TITLE_DICT[plot_name.upper()]}: 1C over time.')
	plt.xlim([0, 360])
	plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/Hart_{plot_name}_1c.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
	plt.cla()
	plt.clf()

	plt.figure(num=1)
	plt.plot(blade_twist_azimuth, cos_pitch_array_hp, linewidth=0.5)
	plt.title(f'{TITLE_DICT[plot_name.upper()]}: 1C over time HP.')
	plt.xlim([0, 360])
	plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/Hart_{plot_name}_1c_hp.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
	plt.cla()
	plt.clf()

	hhc_pitch_array = blade_results['hhc_pitch_array'][2,:]
	
	#print(f'blade_twist_array: {blade_twist_array}')
	plt.figure(num=1)
	plt.plot(blade_twist_azimuth, hhc_pitch_array, linewidth=0.5)
	plt.title(f'{TITLE_DICT[plot_name.upper()]}: HHC over time.')
	plt.xlim([0, 360])
	plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/Hart_{plot_name}_hhc.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
	plt.cla()
	plt.clf()

	blade_flapping_array = blade_results['blade_flapping_array'][2,:]
	
	#print(fblade_twist_array: {blade_twist_array}')
	plt.figure(num=1)
	plt.plot(blade_twist_azimuth, blade_flapping_array*(180.0/math.pi), linewidth=0.5)
	plt.title(f'{TITLE_DICT[plot_name.upper()]}: Flapping over time.')
	plt.xlim([0, 360])
	plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/Hart_{plot_name}_flap.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
	plt.cla()
	plt.clf()

	blade_flapping_der_array = blade_results['blade_flapping_der_array'][2,:]
	
	#print(f'blade_twist_array: {blade_twist_array}')
	plt.figure(num=1)
	plt.plot(blade_twist_azimuth, blade_flapping_der_array, linewidth=0.5)
	plt.title(f'{TITLE_DICT[plot_name.upper()]}: Flapping derivative over time.')
	plt.xlim([0, 360])
	plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/Hart_{plot_name}_flap_dot.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
	plt.cla()
	plt.clf()


if __name__ == "__main__":

	matplotlib.rcParams['font.family'] = 'serif'
	matplotlib.rcParams['font.serif'] = prop.get_name()
	matplotlib.rcParams['mathtext.fontset'] = 'stix'

	plot_blade_twist('BL')
	plot_blade_twist('MN')
	plot_blade_twist('MV')	

	plot_blade_normal_pressures('BL')
	plot_wake_trajectory('BL')
	plot_acoustic_contours("BL")

	plot_blade_normal_pressures('MN')
	plot_wake_trajectory('MN')
	plot_acoustic_contours("MN")

	plot_blade_normal_pressures('MV')
	plot_wake_trajectory('MV')
	plot_acoustic_contours("MV")
