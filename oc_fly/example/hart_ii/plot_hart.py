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

font_size0 = 22
font_size1 = 17
font_size2 = 16
font_size3 = 12
font_size35 = 8
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

FLAPPING_AZIMUTH = {
	"BL": [
		-11.0,
		19.2,
		34.2,
		49.3,
		64.3,
		79.3,
		94.1,
		109.1,
		124.0,
		139.1,
		154.1,
		169.0,
		184.1,
		199.0,
		214.1,
		229.0,
		244.1,
		259.0,
		274.0,
		289.0,
		304.0,
		319.0,
		334.0,
		349.0,
		379.2
	],
	"MN": [
		-10.9,
		19.1,
		34.0,
		49.1,
		64.1,
		79.1,
		94.1,
		109.0,
		124.1,
		139.0,
		154.1,
		169.1,
		184.1,
		199.0,
		214.1,
		229.1,
		244.1,
		259.1,
		274.1,
		289.1,
		304.1,
		319.1,
		334.0,
		349.1,
		379.1
	],
	"MV": [
		-11.0,
		19.1,
		34.1,
		49.1,
		64.0,
		79.1,
		94.0,
		109.1,
		124.0,
		139.0,
		154.0,
		169.0,
		184.0,
		199.0,
		214.0,
		229.0,
		244.1,
		259.0,
		274.1,
		289.0,
		304.1,
		319.0,
		334.0,
		349.0,
		379.1
	]
}

FLAPPING_ANGLE = {
	"BL":
		np.arctan(2/(1.575*100)*
			np.asarray([
				-0.463,
				-0.777,
				-0.836,
				-0.933,
				-1.141,
				-1.356,
				-1.418,
				-1.457,
				-1.337,
				-1.181,
				-0.989,
				-0.827,
				-0.716,
				-0.643,
				-0.466,
				-0.251,
				-0.052,
				0.157,
				0.393,
				0.405,
				0.271,
				0.020,
				-0.260,
				-0.463,
				-0.777
			])
		),
	"MN": np.arctan(2/(1.575*100)*
			np.asarray([
				-0.906,
				-0.683,
				-0.403,
				-0.299,
				-0.584,
				-1.169,
				-1.728,
				-2.095,
				-2.060,
				-1.681,
				-1.227,
				-0.871,
				-0.711,
				-0.666,
				-0.650,
				-0.443,
				-0.031,
				 0.470,
				 0.890,
				 0.927,
				 0.574,
				-0.049,
				-0.627,
				-0.906,
				-0.683
			])
		),
	"MV": np.arctan(2/(1.575*100)*
			np.asarray([
				0.055,
				-0.139,
				-0.614,
				-1.263,
				-1.827,
				-1.993,
				-1.795,
				-1.383,
				-0.910,
				-0.628,
				-0.588,
				-0.691,
				-0.767,
				-0.672,
				-0.434,
				-0.005,
				 0.322,
				 0.400,
				 0.256,
				 0.033,
				-0.181,
				-0.189,
				-0.081,
				 0.055,
				-0.139
			])
		)
}

ELASTIC_TWIST_AZIMUTH = {
	"BL": [
		 19.225,
		 34.245,
		 49.265,
		 64.275,
		 79.285,
		 94.055,
		109.055,
		124.045,
		139.055,
		154.055,
		169.045,
		199.035,
		214.055,
		229.035,
		244.065,
		259.045,
		274.045,
		289.045,
		304.035,
		319.045,
		334.045,
		349.045
	],
	"MN": [
		 19.055,
		 34.045,
		 49.055,
		 64.055,
		 79.055,
		 94.055,
		109.045,
		124.055,
		139.045,
		154.055,
		169.055,
		199.045,
		214.055,
		229.055,
		244.055,
		259.055,
		274.055,
		289.065,
		304.055,
		319.065,
		334.045,
		349.055
	],
	"MV": [
		 19.055,
		 34.055,
		 49.055,
		 64.045,
		 79.055,
		 94.045,
		109.055,
		124.045,
		139.045,
		154.035,
		169.035,
		199.045,
		214.045,
		229.045,
		244.055,
		259.045,
		274.055,
		289.045,
		304.055,
		319.045,
		334.035,
		349.045
	]
}

ELASTIC_TWIST_ANGLE = {
	"BL":
		np.asarray([
			-0.74,
			-0.64,
			-0.63,
			-0.53,
			-0.15,
			-0.01,
			-0.23,
			-0.47,
			-1.06,
			-1.23,
			-1.32,
			-0.93,
			-0.77,
			-0.39,
			-0.11,
			 0.15,
			 0.30,
			 0.24,
			-0.03,
			-0.39,
			-0.66,
			-0.83
		]),
	"MN": np.asarray([
			 0.35,
			-0.50,
			-1.52,
			-1.76,
			-1.14,
			-0.35,
			-0.01,
			 0.55,
			 0.40,
			-0.62,
			-1.72,
			-2.18,
			-1.11,
			 0.11,
			 1.02,
			 1.16,
			 0.48,
			-0.60,
			-1.34,
			-1.33,
			-0.63,
			 0.08
		]),
	"MV": np.asarray([
			-2.24,
			-1.92,
			-1.02,
			 0.21,
			 1.53,
			 1.72,
			 0.73,
			-0.83,
			-2.14,
			-2.70,
			-2.09,
			 0.18,
			 0.68,
			 0.44,
			-0.35,
			-1.03,
			-1.15,
			-0.35,
			 0.63,
			 0.86,
			 0.48,
			-0.52
		]),
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

def read_hart_microphone_data_tecplot(filename: str):
	with open(filename, 'r') as tecfile:
		_ = tecfile.readline()
		_ = tecfile.readline()
		_ = tecfile.readline()
		_ = tecfile.readline()

		#zone_header = tecfile.readline().strip()
		#zone_dictionary = {s[0]: ''.join(s[1:]) for s in [s.split("=") for s in [s.strip() for s in zone_header[4:].split(',')]]}

		num_points = 2048#int(zone_dictionary["I"])
		blade_data = np.zeros((2, num_points))

		for p_idx in range(num_points):
			blade_data[:, p_idx] = [float(x) for x in tecfile.readline().strip().split(' ')]

	return blade_data

cmap_lines = matplotlib.colors.LinearSegmentedColormap.from_list("", ["black", "black"])
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#648FFF", "#785EF0", "#DC267F", "#FE6100", "#FFB000"])
#cmap = matplotlib.colormaps.get_cmap('viridis')

def plot_spectrum(plot_name: str):
	wopwop_results = parse_wopwop_results(f'{os.path.dirname(os.path.realpath(__file__))}/{plot_name.upper()}/acoustics/full_system', 'case.nam')

	mic_data = read_hart_microphone_data_tecplot(f'{os.path.dirname(os.path.realpath(__file__))}/bl_adv-M11-th.tec')

	t_measured = np.linspace(0, 2.0*math.pi/109.12, mic_data.shape[1])
	p_measured = mic_data[1,:]

	bpf = 17.36698741*4
	bvi_spl_start = bpf*6
	bvi_spl_end = bpf*20
	#y = np.asarray(wopwop_results.oaspl_dba_grid.obs_y[0])
	y = np.asarray([_y[0] for _y in wopwop_results.oaspl_dba_grid.obs_y])
	#x = [_x[0] for _x in wopwop_results.oaspl_dba_grid.obs_x]
	x = wopwop_results.oaspl_dba_grid.obs_x[0]
	x_delta = 4 - x[-1]
	x.reverse()
	x = np.asarray(x)
	x = x + x_delta

	i_max = len(wopwop_results.oaspl_db_grid.obs_x)
	j_max = len(wopwop_results.oaspl_db_grid.obs_x[0])

	mic_x = -0.054*2
	#mic_x = 0.1*2
	mic_y = 0.905*2
	#mic_y = 1.1*2

	j = np.argmin(np.abs(x - mic_x))
	i = np.argmin(np.abs(y - mic_y))

	t = np.asarray(wopwop_results.observer_pressures[i*j_max + j].independent_axis)
	p = np.asarray(wopwop_results.observer_pressures[i*j_max + j].functions[2].data)

	t = t - t[0]

	fs = 1.0/(t[1] - t[0])
	fs_measured = 1.0/(t_measured[1] - t_measured[0])

	f, pxx = sig.welch(p, fs, detrend=False, scaling='spectrum', nperseg=len(p))
	f_measured, pxx_measured = sig.welch(p_measured, fs_measured, detrend=False, scaling='spectrum', nperseg=len(p_measured))

	spl = 10*np.log10(pxx/2.0e-5**2)
	spl_measured = 10*np.log10(pxx_measured/2.0e-5**2)
	# p_f = sft.fft(p)
	# p_f_measured = sft.fft(p_measured)
	
	# pxx = 20*np.log10(np.abs(p_f)/2.0e-5)
	# pxx_measured = 20*np.log10(np.abs(p_f_measured)/2.0e-5)

	# f = sft.fftfreq(p.size, t[1] - t[0])
	# f_measured = sft.fftfreq(p_measured.size, t_measured[1] - t_measured[0])

	plt.figure(num = 1)
	plt.plot(f/bpf, pxx, 'b', f_measured/bpf, pxx_measured, 'r.-', linewidth=0.5, markersize=0.7)
	#plt.plot(f_measured, pxx_measured)
	plt.xlim(bvi_spl_start/bpf, bvi_spl_end/bpf)
	#plt.ylim(120, 170)
	#plt.xlim(0, int(bvi_spl_end/2))
	plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/Hart_{plot_name.upper()}_spectrum.pdf', dpi=500, bbox_inches="tight", pad_inches=0.0)
	plt.cla()
	plt.clf()

	plt.figure(num = 1)
	plt.plot(t, p, 'b', t_measured, p_measured, 'r.-', linewidth=0.5, markersize=0.7)
	#plt.xlim(bvi_spl_start, bvi_spl_end)
	plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/Hart_{plot_name.upper()}_acoustic_pressure.pdf', dpi=500, bbox_inches="tight", pad_inches=0.0)
	plt.cla()
	plt.clf()

	#wopwop_results.observer_pressures

def plot_acoustic_contours(plot_name: str):
	x_grid, y_grid, measured = read_hart_contour_tecplot(f'{os.path.dirname(os.path.realpath(__file__))}/{plot_name.lower()}-contour-meas.tec')

	wopwop_results = parse_wopwop_results(f'{os.path.dirname(os.path.realpath(__file__))}/{plot_name.upper()}/acoustics/full_system', 'case.nam')

	i_max = len(wopwop_results.oaspl_db_grid.obs_x)
	j_max = len(wopwop_results.oaspl_db_grid.obs_x[0])

	phi = np.linspace(0, 2*math.pi, 1000)
	rotor_x = 1*np.cos(phi)
	rotor_z = 1*np.sin(phi)

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

	mic_x = -0.054
	#mic_x = 0.1
	mic_y = 0.905
	#mic_y = 1.1


	fig = plt.figure()
	ax0 = plt.subplot(121)
	plt.plot(rotor_x, rotor_z, 'k', linewidth=1.5)
	plt2 = plt.contour([-_y/R for _y in y], [-(_x - offset)/R for _x in x], oaspl_db, levels=clevels, cmap=cmap_lines, linewidths=0.5)
	plt1 = plt.contourf([-_y/R for _y in y], [-(_x - offset)/R for _x in x], oaspl_db, levels=clevels, cmap=cmap)
	#plt.plot(mic_y, mic_x, 'k.', markersize=7)
	#plt.clabel(plt2, clevels, inline=True, colors='k', fontsize=5)
	clabels = plt.clabel(plt2, clevels, colors='k', fontsize=font_size35)
	plt.ylabel('$x/R$', **label_font)
	plt.xlabel('$y/R$', labelpad=20, **label_font)
	plt.title('Prediction', **title_font)
	plt.ylim(-2, 2)
	#plt.axis('equal')
	
	plt.axis('scaled')
	
	for label in clabels:
		label.set_bbox(dict(facecolor='white', edgecolor='white', pad=0.01))
	
	ax = plt.gca()
	for label in (ax.get_xticklabels() + ax.get_yticklabels()):
		label.set_fontname('Times New Roman')
		label.set_fontsize(font_size3)

	ax1 = plt.subplot(122)
	ax1.set_yticklabels([])
	#ax.set_xticklabels([])
	plt.plot(rotor_x, rotor_z, 'k', linewidth=1.5)
	plt2 = plt.contour(x_grid[0,:], y_grid[:,0], measured, levels=clevels, cmap=cmap_lines, linewidths=0.5)
	plt.contourf(x_grid[0,:], y_grid[:,0], measured, levels=clevels, cmap=cmap)
	#plt.plot(mic_y, mic_x, 'k.', markersize=7)
	clabels = plt.clabel(plt2, clevels, colors='k', fontsize=font_size35)
	plt.title('Measured', **title_font)
	#plt.clabel(CS, clevels, inline=True)
	#plt.contourf(y_grid[:,0], x_grid[0,:], measured, levels=clevels)
	
	#plt.axis('equal')

	plt.ylim(-2, 2)
	plt.axis('scaled')
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

	plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/Hart_{plot_name.upper()}.pdf', dpi=500, bbox_inches="tight", pad_inches=0.0)
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
	#plt.title(f'{TITLE_DICT[plot_name.upper()]} 87% Span blade loading. Mean removed.')
	plt.legend(['Measured', 'OpenCOPTER'], ncol=2, frameon=False, fancybox=False)
	plt.xlabel("$\psi$")
	plt.ylabel("$C_nM^2$")
	plt.ylim(-0.2, 0.2)
	plt.xticks(ticks=np.linspace(0, 360, 5))
	plt.gca().set_aspect(400)
	plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/Hart_{plot_name}_normal_loading.pdf', dpi=500, bbox_inches="tight", pad_inches=0.01)
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
	#plt.title(f'{TITLE_DICT[plot_name.upper()]}: 87% Span blade loading. 10/rev highpassed.')
	plt.legend(['Measured', 'OpenCOPTER'], ncol=2, frameon=False, fancybox=False)
	#plt.legend(['Measured', 'OpenCOPTER'])
	#plt.xlim([20, 90])
	plt.ylim(-0.05, 0.05)
	plt.xticks(ticks=np.linspace(0, 360, 5))
	plt.gca().set_aspect(1000)
	plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/Hart_{plot_name}_normal_loading_hp.pdf', dpi=500, bbox_inches="tight", pad_inches=0.01)
	plt.cla()
	plt.clf()

	measured_adv_side_range = np.logical_and(measured_azimuth < 90, measured_azimuth > 10)
	comp_adv_side_range = np.logical_and(computational_azimuth < 90, computational_azimuth > 10)

	measured_ret_side_range = np.logical_and(measured_azimuth < 360, measured_azimuth > 280)
	comp_ret_side_range = np.logical_and(computational_azimuth < 360, computational_azimuth > 280)
	
	plt.figure(num=1)
	plt.plot(measured_azimuth[measured_adv_side_range], measured_blade_loading[measured_adv_side_range], 'r.-', computational_azimuth[comp_adv_side_range], blade_loading[comp_adv_side_range], 'b', linewidth=0.5, markersize=0.7)
	plt.xlabel("$\psi$", fontsize=font_size0)
	plt.ylabel("$C_nM^2$", fontsize=font_size0)
	#plt.title(f'{TITLE_DICT[plot_name.upper()]}: 87% Span blade loading. 10/rev highpassed. Advancing side')
	plt.legend(['Measured', 'OpenCOPTER'], ncol=2, frameon=False, fancybox=False, fontsize=font_size1)
	#plt.legend(['Measured', 'OpenCOPTER'])
	#plt.xlim([20, 90])
	plt.ylim(-0.05, 0.05)
	plt.xticks(ticks=np.linspace(10, 90, 9), fontsize=font_size0)
	plt.yticks(fontsize=font_size0)
	plt.gca().set_aspect(300)
	plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/Hart_{plot_name}_normal_loading_hp_adv.pdf', dpi=500, bbox_inches="tight", pad_inches=0.01)
	plt.cla()
	plt.clf()

	plt.figure(num=1)
	plt.plot(measured_azimuth[measured_ret_side_range], measured_blade_loading[measured_ret_side_range], 'r.-', computational_azimuth[comp_ret_side_range], blade_loading[comp_ret_side_range], 'b', linewidth=0.5, markersize=0.7)
	plt.xlabel("$\psi$", fontsize=font_size0)
	#plt.ylabel("$C_nM^2$", fontsize=font_size0)
	#plt.title(f'{TITLE_DICT[plot_name.upper()]}: 87% Span blade loading. 10/rev highpassed. Retreating side')
	plt.legend(['Measured', 'OpenCOPTER'], ncol=2, frameon=False, fancybox=False, fontsize=font_size1)
	#plt.legend(['Measured', 'OpenCOPTER'])
	#plt.xlim([20, 90])
	plt.ylim(-0.05, 0.05)
	plt.xticks(ticks=np.linspace(280, 360, 9), fontsize=font_size0)
	plt.yticks(ticks=[], fontsize=font_size0)
	plt.gca().set_aspect(300)
	plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/Hart_{plot_name}_normal_loading_hp_ret.pdf', dpi=500, bbox_inches="tight", pad_inches=0.01)
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
		plt.legend(['Measured', 'OpenCOPTER'], ncol=2, frameon=False, fancybox=False)
		#plt.legend(['Measured', 'OpenCOPTER'])
		#plt.axis('scaled')
		plt.xlim([-1.1, 1.1])
		plt.ylim([-0.1, 0.15])
		plt.gca().set_aspect(3)
		plt.xlabel("$x/R$")
		plt.ylabel("$z/R$")
		#plt.gca().invert_xaxis()
		plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/Hart_{plot_name}_wake_y{target_y_slice:.2f}.pdf', dpi=500, bbox_inches="tight", pad_inches=0.01)
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
		plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/Hart_{plot_name}_core_size_y{target_y_slice:.2f}.pdf', dpi=500, bbox_inches="tight", pad_inches=0.01)
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
	plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/Hart_{plot_name}_fil_core_size.pdf', dpi=500, bbox_inches="tight", pad_inches=0.0)
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

	#blade_twist_array = blade_results['blade_twist_array'][2,:]
	blade_twist_azimuth = blade_results['blade_twist_azimuth'][2,:]
	

	# computed_dt = (2.0*math.pi/omega)/blade_twist_array.shape[0]
	# computed_sos = sig.butter(13, ten_per_rev, 'highpass', output='sos', fs=1/computed_dt)

	# #print(f'blade_twist_array: {blade_twist_array}')
	# plt.figure(num=1)
	# plt.plot(blade_twist_azimuth, blade_twist_array, linewidth=0.5)
	# plt.title(f'{TITLE_DICT[plot_name.upper()]}: Blade twist over time.')
	# plt.xlim([0, 360])
	# plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/Hart_{plot_name}_blade_twist.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
	# plt.cla()
	# plt.clf()
	
	elastic_twist_array = blade_results['elastic_twist_array'][2,:]
	
	#print(f'blade_twist_array: {blade_twist_array}')
	plt.figure(num=1)
	plt.plot(blade_twist_azimuth, elastic_twist_array*(180.0/math.pi), 'b', linewidth=1.0)
	plt.plot(ELASTIC_TWIST_AZIMUTH[plot_name], ELASTIC_TWIST_ANGLE[plot_name], 'r*', linewidth=1.0, markersize=10)
	#plt.title(f'{TITLE_DICT[plot_name.upper()]}: Elastic twist over time.')
	plt.xlabel('$\psi$ [$^\circ$]', fontsize=font_size0)
	
	if plot_name == "BL":
		plt.ylabel('$\\theta_e$ [$^\circ$]', fontsize=font_size0)
		plt.legend(['Fit', 'Measured'], ncol=2, frameon=False, fancybox=False, fontsize=font_size0)
		plt.yticks(fontsize=font_size0)
	else:
		plt.yticks(ticks = [], fontsize=font_size0)

	plt.xlim([0, 360])
	plt.ylim([-3, 2])
	plt.gca().set_aspect(60)
	plt.xticks(ticks=np.linspace(0, 360, 5), fontsize=font_size0)
	
	plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/Hart_{plot_name}_elastic_twist.pdf', dpi=500, bbox_inches="tight", pad_inches=0.01)
	plt.cla()
	plt.clf()

	# collective_pitch_array = blade_results['collective_pitch_array'][:]
	
	# #print(f'blade_twist_array: {blade_twist_array}')
	# plt.figure(num=1)
	# plt.plot(blade_twist_azimuth, collective_pitch_array, linewidth=0.5)
	# plt.title(f'{TITLE_DICT[plot_name.upper()]}: Collective pitch over time.')
	# plt.xlim([0, 360])
	# plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/Hart_{plot_name}_collective.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
	# plt.cla()
	# plt.clf()

	# sin_pitch_array = blade_results['sin_pitch_array'][2,:]
	# sin_pitch_array_hp = sig.sosfilt(computed_sos, sin_pitch_array)

	# #print(f'blade_twist_array: {blade_twist_array}')
	# plt.figure(num=1)
	# plt.plot(blade_twist_azimuth, sin_pitch_array, linewidth=0.5)
	# plt.title(f'{TITLE_DICT[plot_name.upper()]}: 1S over time.')
	# plt.xlim([0, 360])
	# plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/Hart_{plot_name}_1s.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
	# plt.cla()
	# plt.clf()

	# #print(f'blade_twist_array: {blade_twist_array}')
	# plt.figure(num=1)
	# plt.plot(blade_twist_azimuth, sin_pitch_array_hp, linewidth=0.5)
	# plt.title(f'{TITLE_DICT[plot_name.upper()]}: 1S over time HP.')
	# plt.xlim([0, 360])
	# plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/Hart_{plot_name}_1s_hp.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
	# plt.cla()
	# plt.clf()

	# cos_pitch_array = blade_results['cos_pitch_array'][2,:]
	# cos_pitch_array_hp = sig.sosfilt(computed_sos, cos_pitch_array)

	# #print(f'blade_twist_array: {blade_twist_array}')
	# plt.figure(num=1)
	# plt.plot(blade_twist_azimuth, cos_pitch_array, linewidth=0.5)
	# plt.title(f'{TITLE_DICT[plot_name.upper()]}: 1C over time.')
	# plt.xlim([0, 360])
	# plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/Hart_{plot_name}_1c.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
	# plt.cla()
	# plt.clf()

	# plt.figure(num=1)
	# plt.plot(blade_twist_azimuth, cos_pitch_array_hp, linewidth=0.5)
	# plt.title(f'{TITLE_DICT[plot_name.upper()]}: 1C over time HP.')
	# plt.xlim([0, 360])
	# plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/Hart_{plot_name}_1c_hp.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
	# plt.cla()
	# plt.clf()

	# hhc_pitch_array = blade_results['hhc_pitch_array'][2,:]
	
	# #print(f'blade_twist_array: {blade_twist_array}')
	# plt.figure(num=1)
	# plt.plot(blade_twist_azimuth, hhc_pitch_array, linewidth=0.5)
	# plt.title(f'{TITLE_DICT[plot_name.upper()]}: HHC over time.')
	# plt.xlim([0, 360])
	# plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/Hart_{plot_name}_hhc.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
	# plt.cla()
	# plt.clf()

	blade_flapping_array = blade_results['blade_flapping_array'][2,:]
	
	#print(fblade_twist_array: {blade_twist_array}')
	plt.figure(num=1)
	plt.plot(blade_twist_azimuth, blade_flapping_array*(180.0/math.pi), 'b', linewidth=1.0)
	plt.plot(FLAPPING_AZIMUTH[plot_name], FLAPPING_ANGLE[plot_name]*(180.0/math.pi), 'r*', linewidth=1.0, markersize=10)
	#plt.title(f'{TITLE_DICT[plot_name.upper()]}: Flapping over time.')
	#if plot_name == "BL":
	#	plt.legend(['Fit', 'Measured'], ncol=2, frameon=False, fancybox=False)
	if plot_name == "BL":
		plt.ylabel('$\\beta$ [$^\circ$]', fontsize=font_size0)
		plt.legend(['Fit', 'Measured'], ncol=2, frameon=False, fancybox=False, fontsize=font_size0)
		plt.yticks(fontsize=font_size0)
	else:
		plt.yticks(ticks = [], fontsize=font_size0)
	plt.xlim([0, 360])
	plt.ylim([-1.6, 0.8])
	plt.xlabel('$\psi$ [$^\circ$]', fontsize=font_size0)
	#plt.ylabel('$\\beta$ [$^\circ$]', fontsize=font_size0)
	plt.xticks(ticks=np.linspace(0, 360, 5), fontsize=font_size0)
	#plt.yticks(fontsize=font_size0)
	plt.gca().set_aspect(120)
	plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/Hart_{plot_name}_flap.pdf', dpi=500, bbox_inches="tight", pad_inches=0.01)
	plt.cla()
	plt.clf()

	# blade_flapping_der_array = blade_results['blade_flapping_der_array'][2,:]
	
	# #print(f'blade_twist_array: {blade_twist_array}')
	# plt.figure(num=1)
	# plt.plot(blade_twist_azimuth, blade_flapping_der_array, linewidth=0.5)
	# plt.title(f'{TITLE_DICT[plot_name.upper()]}: Flapping derivative over time.')
	# plt.xlim([0, 360])
	# plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/Hart_{plot_name}_flap_dot.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
	# plt.cla()
	# plt.clf()


if __name__ == "__main__":

	matplotlib.rcParams['font.family'] = 'serif'
	matplotlib.rcParams['font.serif'] = prop.get_name()
	matplotlib.rcParams['mathtext.fontset'] = 'stix'

	#plot_spectrum('BL')

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
