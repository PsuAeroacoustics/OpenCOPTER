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
import csv

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
font_size4 = 7
font_size5 = 5

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
	"ID1": "ID1",
	"ID1_ASD": "ID1 Advancing side down",
	"ID1_ASU": "ID1 Advancing side up",
	"ID2": "ID2",
	"ID2_ASD": "ID2 Advancing side down",
	"ID2_ASU": "ID2 Advancing side up",
	"ID5": "ID5"
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

def read_hart_blade_data_tsv(filename: str):
	azimuth = []
	loading = []
	with open(filename, 'r') as tsv_file:

		loading_reader = csv.reader(tsv_file, delimiter='\t')
		for row in loading_reader:
			azimuth.append(float(row[0]))
			loading.append(float(row[1]))

	blade_data = np.zeros((2, len(azimuth)))
	blade_data[0,:] = azimuth
	blade_data[1,:] = loading

	return blade_data

ACOUSTIC_BOUNDS_DB = {
	"ID1": [(82, 100, 19), (82, 100, 19), (81, 109, 29)],
	"ID1_ASU": [(82, 100, 19), (82, 100, 19), (81, 109, 29)],
	"ID1_ASD": [(82, 100, 19), (82, 100, 19), (81, 109, 29)],
	"ID2": [(81, 99, 19), (81, 99, 19), (81, 109, 29)],
	"ID2_ASU": [(81, 99, 19), (81, 99, 19), (81, 109, 29)],
	"ID2_ASD": [(81, 99, 19), (81, 99, 19), (81, 109, 29)],
	"ID5": [(85, 119, 18), (85, 119, 18), (85, 119, 18)]
}

ACOUSTIC_BOUNDS_DBA = {
	"ID1": [(82, 100, 19), (82, 100, 19), (84, 101, 18)],
	"ID1_ASU": [(82, 100, 19), (82, 100, 19), (82, 100, 19)],
	"ID1_ASD": [(82, 100, 19), (82, 100, 19), (82, 100, 19)],
	"ID2": [(81, 99, 19), (81, 99, 19), (84, 104, 21)],
	"ID2_ASU": [(81, 99, 19), (81, 99, 19), (81, 99, 19)],
	"ID2_ASD": [(81, 99, 19), (81, 99, 19), (81, 99, 19)],
	"ID5": [(85, 119, 18), (85, 119, 18), (96, 113, 18)]
}

ACOUSTIC_BOUNDS_TR_DB = {
	"ID1": [(82, 100, 19), (82, 100, 19), (81, 109, 29)],
	"ID1_ASU": [(82, 100, 19), (82, 100, 19), (81, 109, 29)],
	"ID1_ASD": [(82, 100, 19), (82, 100, 19), (81, 109, 29)],
	"ID2": [(81, 99, 19), (81, 99, 19), (81, 109, 29)],
	"ID2_ASU": [(81, 99, 19), (81, 99, 19), (81, 109, 29)],
	"ID2_ASD": [(81, 99, 19), (81, 99, 19), (81, 109, 29)],
	"ID5": [(85, 119, 18), (85, 119, 18), (85, 119, 18)]
}

ACOUSTIC_BOUNDS_TR_DBA = {
	"ID1": [(82, 100, 19), (82, 100, 19), (84, 101, 18)],
	"ID1_ASU": [(82, 100, 19), (82, 100, 19), (84, 101, 18)],
	"ID1_ASD": [(82, 100, 19), (82, 100, 19), (84, 101, 18)],
	"ID2": [(81, 99, 19), (81, 99, 19), (81, 99, 19)],
	"ID2_ASU": [(81, 99, 19), (81, 99, 19), (81, 99, 19)],
	"ID2_ASD": [(81, 99, 19), (81, 99, 19), (81, 99, 19)],
	"ID5": [(85, 119, 18), (85, 119, 18), (85, 119, 18)]
}

#cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["b", "dodgerblue", "cyan", "chartreuse", "yellow", "gold", "orange", "darkorange", "red"])
#cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["b", "dodgerblue", "cyan", "springgreen", "chartreuse", "greenyellow","yellow", "orange", "red"])


#cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["b", "dodgerblue", "cyan", "springgreen", "chartreuse", "yellow", "orange", "red"])
#cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["b", "dodgerblue", "cyan", "springgreen", "chartreuse", "yellow", "orange", "red", "crimson", "magenta"])#, "purple"])
#cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["b", "dodgerblue", "cyan", "springgreen", "chartreuse", "yellow", "orange", "red", "magenta"])#, "purple"])
cmap_lines = matplotlib.colors.LinearSegmentedColormap.from_list("", ["black", "black"])

ACOUSTIC_CMAP = {
	"ID1": matplotlib.colors.LinearSegmentedColormap.from_list("", ["b", "dodgerblue", "cyan", "springgreen", "chartreuse", "yellow", "orange", "red"]),
	"ID1_ASU": matplotlib.colors.LinearSegmentedColormap.from_list("", ["b", "dodgerblue", "cyan", "springgreen", "chartreuse", "yellow", "orange", "red"]),
	"ID1_ASD": matplotlib.colors.LinearSegmentedColormap.from_list("", ["b", "dodgerblue", "cyan", "springgreen", "chartreuse", "yellow", "orange", "red"]),
	"ID2": matplotlib.colors.LinearSegmentedColormap.from_list("", ["b", "dodgerblue", "cyan", "springgreen", "chartreuse", "yellow", "orange", "red", "magenta"]),
	"ID2_ASU": matplotlib.colors.LinearSegmentedColormap.from_list("", ["b", "dodgerblue", "cyan", "springgreen", "chartreuse", "yellow", "orange", "red"]),
	"ID2_ASD": matplotlib.colors.LinearSegmentedColormap.from_list("", ["b", "dodgerblue", "cyan", "springgreen", "chartreuse", "yellow", "orange", "red"]),
	"ID5": matplotlib.colors.LinearSegmentedColormap.from_list("", ["b", "dodgerblue", "cyan", "springgreen", "chartreuse", "yellow", "orange", "red"]),
}

MEASURED_DATA_FS_DBA = {
	"ID1": ["", "", "HeliNOVI_ID1_Full_spectrum_dBA_measured.jpg"],
	"ID1_ASU": ["", "", ""],
	"ID1_ASD": ["", "", ""],
	"ID2": ["", "", "HeliNOVI_ID2_Full_spectrum_dBA_measured.jpg"],
	"ID2_ASU": ["", "", ""],
	"ID2_ASD": ["", "", ""],
	"ID5": ["", "", "HeliNOVI_ID5_Full_spectrum_dBA_measured.jpg"],
}

UPM_DATA_FS_DBA = {
	"ID1": ["", "", "HeliNOVI_ID1_Full_spectrum_dBA_upm.jpg"],
	"ID1_ASU": ["", "", ""],
	"ID1_ASD": ["", "", ""],
	"ID2": ["", "", "HeliNOVI_ID2_Full_spectrum_dBA_upm.jpg"],
	"ID2_ASU": ["", "", ""],
	"ID2_ASD": ["", "", ""],
	"ID5": ["", "", "HeliNOVI_ID5_Full_spectrum_dBA_upm.jpg"],
}

MEASURED_DATA_TR_DBA = {
	"ID1": ["", "", ""],
	"ID1_ASU": ["", "", "HeliNOVI_ID1_ASU_Full_spectrum_TR_dBA_measured.jpg"],
	"ID1_ASD": ["", "", "HeliNOVI_ID1_ASD_Full_spectrum_TR_dBA_measured.jpg"],
	"ID2": ["", "", ""],
	"ID2_ASU": ["", "", ""],
	"ID2_ASD": ["", "", ""],
	"ID5": ["", "", ""],
}

UPM_DATA_TR_DBA = {
	"ID1": ["", "", ""],
	"ID1_ASU": ["", "", "HeliNOVI_ID1_ASU_Full_spectrum_TR_dBA_upm.jpg"],
	"ID1_ASD": ["", "", "HeliNOVI_ID1_ASD_Full_spectrum_TR_dBA_upm.jpg"],
	"ID2": ["", "", ""],
	"ID2_ASU": ["", "", ""],
	"ID2_ASD": ["", "", ""],
	"ID5": ["", "", ""],
}

def plot_acoustic_contours_fs(plot_name: str, plot_hybrid = False):
	#x_grid, y_grid, measured = read_hart_contour_tecplot(f'{os.path.dirname(os.path.realpath(__file__))}/{plot_name.lower()}-contour-meas.tec')

	namelist = parse_namelist(f'{os.path.dirname(os.path.realpath(__file__))}/{plot_name.upper()}/acoustics/full_system/case.nam')
	wopwop_results = parse_wopwop_results_namelist(f'{os.path.dirname(os.path.realpath(__file__))}/{plot_name.upper()}/acoustics/full_system', namelist)

	hybrid_wopwop_results = None
	if plot_hybrid:
		hybrid_wopwop_results = parse_wopwop_results_namelist(f'{os.path.dirname(os.path.realpath(__file__))}/../helinovi_hybrid/helinovi/{plot_name.upper()}/acoustics/full_system', namelist)

	i_max = len(wopwop_results.oaspl_dba_grid.obs_x)
	j_max = len(wopwop_results.oaspl_dba_grid.obs_x[0])

	#print(f'i_max: {i_max}')
	#print(f'j_max: {j_max}')

	phi = np.linspace(0, 2*math.pi, 1000)
	rotor_x = 2*np.cos(phi)
	rotor_z = 2*np.sin(phi)

	for r_idx, freq_range in enumerate(namelist.observers[0].ranges):
		oaspl_linear = [oaspl_db.functions[2].data[0] for oaspl_db in wopwop_results.frequency_ranges_db[r_idx]]

		oaspl_db = [[oaspl_linear[j*i_max + i] for i in range(i_max)] for j in range(j_max)]

		oaspl_linear = [oaspl_db.functions[2].data[0] for oaspl_db in wopwop_results.frequency_ranges_dba[r_idx]]

		oaspl_dba = [[oaspl_linear[j*i_max + i] for i in range(i_max)] for j in range(j_max)]

		#print(wopwop_results.oaspl_db_grid.obs_y[0])
		#print(wopwop_results.oaspl_db_grid.obs_x[1])

		y = wopwop_results.oaspl_dba_grid.obs_y[0]
		x = [_x[0] for _x in wopwop_results.oaspl_dba_grid.obs_x]
		x_delta = 4 - x[0]
		x.reverse()
		x = np.asarray(x)
		x = x + x_delta

		clevels = np.linspace(ACOUSTIC_BOUNDS_DB[plot_name][r_idx][0], ACOUSTIC_BOUNDS_DB[plot_name][r_idx][1], ACOUSTIC_BOUNDS_DB[plot_name][r_idx][2])
		cnorm = plt.Normalize(ACOUSTIC_BOUNDS_DB[plot_name][r_idx][0], ACOUSTIC_BOUNDS_DB[plot_name][r_idx][1])

		fig = plt.figure()

		plt2 = plt.contour(y, x, oaspl_db, levels=clevels, linewidths=0.5, cmap=cmap_lines, norm=cnorm)
		plt1 = plt.contourf(y, x, oaspl_db, levels=clevels, cmap=ACOUSTIC_CMAP[plot_name], norm=cnorm)
		clabels = plt.clabel(plt2, clevels, colors='k', fontsize=font_size4)
		plt.ylabel('x [m]')
		plt.xlabel('y [m]')
		plt.title('OpenCOPTER Standard')
		plt.axis('scaled')

		for label in clabels:
			label.set_bbox(dict(facecolor='white', edgecolor='white', pad=0.01))

		cax = fig.add_axes([0.95, 0.11, 0.02, 0.77])
		fig.colorbar(plt1, cax, orientation='vertical', label='BVI SPL [dB]')
		plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name.upper()}_{freq_range.Title}_dB.pdf', dpi=500, bbox_inches="tight", pad_inches=0.0)

		clevels = np.linspace(ACOUSTIC_BOUNDS_DBA[plot_name][r_idx][0], ACOUSTIC_BOUNDS_DBA[plot_name][r_idx][1], ACOUSTIC_BOUNDS_DBA[plot_name][r_idx][2])
		cnorm = plt.Normalize(ACOUSTIC_BOUNDS_DBA[plot_name][r_idx][0], ACOUSTIC_BOUNDS_DBA[plot_name][r_idx][1])

		num_subplots = 1
		if MEASURED_DATA_FS_DBA[plot_name][r_idx] != "":
			num_subplots = num_subplots + 1

		if (UPM_DATA_FS_DBA[plot_name][r_idx] != "") or plot_hybrid:
			num_subplots = num_subplots + 1

		fig = plt.figure()
		_ = plt.subplot(int("1"+str(num_subplots)+"1"))

		plt2 = plt.contour(y, x, oaspl_dba, levels=clevels, linewidths=0.5, cmap=cmap_lines, norm=cnorm)
		plt1 = plt.contourf(y, x, oaspl_dba, levels=clevels, cmap=ACOUSTIC_CMAP[plot_name], norm=cnorm)
		plt.plot(rotor_x, rotor_z, "k", linewidth=1)

		if plot_name != "ID5":
			clabels = plt.clabel(plt2, clevels[0::2], colors='k', fontsize=font_size4)
		else:
			clabels = plt.clabel(plt2, clevels[1::2], colors='k', fontsize=font_size4)
		plt.ylabel('x [m]')
		plt.xlabel('y [m]')
		
		if plot_hybrid:
			plt.title('OpenCOPTER Standard')
		else:
			plt.title('Prediction: OpenCOPTER')
			
		plt.axis('scaled')
		plt.ylim(-4, 4)

		for label in clabels:
			label.set_bbox(dict(facecolor='white', edgecolor='white', pad=0.01))

		plotted_subplots = 2
		if MEASURED_DATA_FS_DBA[plot_name][r_idx] != "":
			_ = plt.subplot(int("1"+str(num_subplots)+str(plotted_subplots)))
			arr = plt.imread(f'{os.path.dirname(os.path.realpath(__file__))}/{MEASURED_DATA_FS_DBA[plot_name][r_idx]}')

			plt.imshow(np.flipud(arr) ,interpolation='bilinear', origin='lower', extent=[y[0],y[-1],-4,4])
			plt.plot(rotor_x, rotor_z, "k", linewidth=1)
			plt.title('Measured')
			plt.xlabel('y [m]')
		
			plotted_subplots = plotted_subplots+1

		if (UPM_DATA_FS_DBA[plot_name][r_idx] != "") and not plot_hybrid:
			_ = plt.subplot(int("1"+str(num_subplots)+str(plotted_subplots)))
			arr = plt.imread(f'{os.path.dirname(os.path.realpath(__file__))}/{UPM_DATA_FS_DBA[plot_name][r_idx]}')
			
			plt.imshow(np.flipud(arr) ,interpolation='bilinear', origin='lower', extent=[y[0],y[-1],-4,4])
			plt.plot(rotor_x, rotor_z, "k", linewidth=1)
			plt.title('Prediction: UPM')
			plt.xlabel('y [m]')
		elif plot_hybrid:
			hybrid_oaspl_linear = [oaspl_db.functions[2].data[0] for oaspl_db in hybrid_wopwop_results.frequency_ranges_dba[r_idx]]

			hybrid_oaspl_dba = [[hybrid_oaspl_linear[j*i_max + i] for i in range(i_max)] for j in range(j_max)]

			_ = plt.subplot(int("1"+str(num_subplots)+str(plotted_subplots)))
			plt2 = plt.contour(y, x, hybrid_oaspl_dba, levels=clevels, linewidths=0.5, cmap=cmap_lines, norm=cnorm)
			plt1 = plt.contourf(y, x, hybrid_oaspl_dba, levels=clevels, cmap=ACOUSTIC_CMAP[plot_name], norm=cnorm)
			plt.plot(rotor_x, rotor_z, "k", linewidth=1)

			clabels = plt.clabel(plt2, clevels[0::2], colors='k', fontsize=font_size4)
			for label in clabels:
				label.set_bbox(dict(facecolor='white', edgecolor='white', pad=0.01))

			plt.xlabel('y [m]')
			
			plt.title('OpenCOPTER Hybrid')
			plt.axis('scaled')
			plt.ylim(-4, 4)
		
		cax = fig.add_axes([0.95, 0.31, 0.02, 0.37])
		fig.colorbar(plt1, cax, orientation='vertical', label='OASPL [dBA]')
		if not plot_hybrid:
			plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name.upper()}_{freq_range.Title}_dBA.pdf', dpi=500, bbox_inches="tight", pad_inches=0.0)
		else:
			plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name.upper()}_{freq_range.Title}_hybrid_dBA.pdf', dpi=500, bbox_inches="tight", pad_inches=0.0)

		plt.clf()
		plt.cla()


def plot_acoustic_contours_mr(plot_name: str):
	#x_grid, y_grid, measured = read_hart_contour_tecplot(f'{os.path.dirname(os.path.realpath(__file__))}/{plot_name.lower()}-contour-meas.tec')

	namelist = parse_namelist(f'{os.path.dirname(os.path.realpath(__file__))}/{plot_name.upper()}/acoustics/rotor_0/case.nam')
	wopwop_results = parse_wopwop_results_namelist(f'{os.path.dirname(os.path.realpath(__file__))}/{plot_name.upper()}/acoustics/rotor_0', namelist)

	i_max = len(wopwop_results.oaspl_dba_grid.obs_x)
	j_max = len(wopwop_results.oaspl_dba_grid.obs_x[0])

	#print(f'i_max: {i_max}')
	#print(f'j_max: {j_max}')

	for r_idx, freq_range in enumerate(namelist.observers[0].ranges):
		oaspl_linear = [oaspl_db.functions[0].data[0] for oaspl_db in wopwop_results.frequency_ranges_db[r_idx]]

		oaspl_db = [[oaspl_linear[j*i_max + i] for i in range(i_max)] for j in range(j_max)]

		oaspl_linear = [oaspl_db.functions[0].data[0] for oaspl_db in wopwop_results.frequency_ranges_dba[r_idx]]

		oaspl_dba = [[oaspl_linear[j*i_max + i] for i in range(i_max)] for j in range(j_max)]

		#print(wopwop_results.oaspl_db_grid.obs_y[0])
		#print(wopwop_results.oaspl_db_grid.obs_x[1])

		y = wopwop_results.oaspl_dba_grid.obs_y[0]
		x = [_x[0] for _x in wopwop_results.oaspl_dba_grid.obs_x]
		x.reverse()

		#print(f'x: {x}')
		#print(f'y: {y}')
		#print(x)
		offset = x[0] - -4.0

		#clevels = np.linspace(85, 119, 18)
		clevels = np.linspace(ACOUSTIC_BOUNDS_DB[plot_name][r_idx][0], ACOUSTIC_BOUNDS_DB[plot_name][r_idx][1], ACOUSTIC_BOUNDS_DB[plot_name][r_idx][2])
		cnorm = plt.Normalize(ACOUSTIC_BOUNDS_DB[plot_name][r_idx][0], ACOUSTIC_BOUNDS_DB[plot_name][r_idx][1])
		#light_rainbow = cmap_map(lambda x: x/2 + 0.5, matplotlib.cm.rainbow)

		#print([(_x - offset)/R for _x in x])
		#print([_y/R for _y in y])
		#print(x_grid[0,:])
		#print(y_grid[:,0])
		#print(x)

		#print(f'oaspl_db: {oaspl_db}')
		fig = plt.figure()
		#ax0 = plt.subplot(111)
		#plt2 = plt.contour([-_y/R for _y in y], [-(_x - offset)/R for _x in x], oaspl_db, levels=clevels, linewidths=0.1)
		#plt1 = plt.contourf([-_y/R for _y in y], [-(_x - offset)/R for _x in x], oaspl_db, levels=clevels)
		#plt2 = plt.contour([-_y/R for _y in y], [-(_x - offset)/R for _x in x], oaspl_db, linewidths=0.1)
		#plt1 = plt.contourf([-_y/R for _y in y], [-(_x - offset)/R for _x in x], oaspl_db)
		#plt2 = plt.contour(y, x, oaspl_db, levels=clevels, linewidths=0.1)
		#plt1 = plt.contourf(y, x, oaspl_db, levels=clevels)
		plt2 = plt.contour(y, x, oaspl_db, levels=15, linewidths=0.5, cmap=cmap_lines, norm=cnorm)
		plt1 = plt.contourf(y, x, oaspl_db, levels=15, cmap=ACOUSTIC_CMAP[plot_name], norm=cnorm)
		#plt.clabel(plt2, clevels, inline=True, colors='k', fontsize=5)
		#plt.clabel(plt2, clevels, colors='k', fontsize=5)
		#plt.colorbar()
		plt.ylabel('x/R')
		plt.xlabel('y/R')
		plt.title('Prediction')
		plt.axis('scaled')
		#plt.axis('equal')

		# ax1 = plt.subplot(122)
		# ax1.set_yticklabels([])
		# #ax.set_xticklabels([])
		# plt2 = plt.contour(x_grid[0,:], y_grid[:,0], measured, levels=clevels, linewidths=0.1)
		# plt.contourf(x_grid[0,:], y_grid[:,0], measured, levels=clevels)
		# plt.clabel(plt2, clevels, colors='k', fontsize=5)
		# plt.title('Measured')
		# #plt.clabel(CS, clevels, inline=True)
		# #plt.contourf(y_grid[:,0], x_grid[0,:], measured, levels=clevels)
		# #plt.axis('equal')
		# #plt.colorbar()
		# plt.xlabel('y/R')
		#plt.ylabel('-x/R')

		#ax = plt.subplot(133)
		cax = fig.add_axes([0.95, 0.11, 0.02, 0.77])
		fig.colorbar(plt1, cax, orientation='vertical', label='BVI SPL [dB]')
		plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name.upper()}_{freq_range.Title}_MR_dB.pdf', dpi=500, bbox_inches="tight", pad_inches=0.0)
		#plt.show()


		clevels = np.linspace(ACOUSTIC_BOUNDS_DBA[plot_name][r_idx][0], ACOUSTIC_BOUNDS_DBA[plot_name][r_idx][1], ACOUSTIC_BOUNDS_DBA[plot_name][r_idx][2])
		cnorm = plt.Normalize(ACOUSTIC_BOUNDS_DBA[plot_name][r_idx][0], ACOUSTIC_BOUNDS_DBA[plot_name][r_idx][1])
		#light_rainbow = cmap_map(lambda x: x/2 + 0.5, matplotlib.cm.rainbow)

		#print([(_x - offset)/R for _x in x])
		#print([_y/R for _y in y])
		#print(x_grid[0,:])
		#print(y_grid[:,0])
		#print(x)

		#print(f'oaspl_db: {oaspl_db}')
		fig = plt.figure()
		#ax0 = plt.subplot(111)
		#plt2 = plt.contour([-_y/R for _y in y], [-(_x - offset)/R for _x in x], oaspl_db, levels=clevels, linewidths=0.1)
		#plt1 = plt.contourf([-_y/R for _y in y], [-(_x - offset)/R for _x in x], oaspl_db, levels=clevels)
		#plt2 = plt.contour([-_y/R for _y in y], [-(_x - offset)/R for _x in x], oaspl_db, linewidths=0.1)
		#plt1 = plt.contourf([-_y/R for _y in y], [-(_x - offset)/R for _x in x], oaspl_db)
		#plt2 = plt.contour(y, x, oaspl_dba, levels=clevels, linewidths=0.1)
		#plt1 = plt.contourf(y, x, oaspl_dba, levels=clevels)
		plt2 = plt.contour(y, x, oaspl_dba, levels=15, linewidths=0.1, cmap=cmap_lines, norm=cnorm)
		plt1 = plt.contourf(y, x, oaspl_dba, levels=15, cmap=ACOUSTIC_CMAP[plot_name], norm=cnorm)
		#plt.clabel(plt2, clevels, inline=True, colors='k', fontsize=5)
		#plt.clabel(plt2, clevels, colors='k', fontsize=5)
		#plt.colorbar()
		plt.ylabel('x/R')
		plt.xlabel('y/R')
		plt.title('Prediction')
		plt.axis('scaled')
		#plt.axis('equal')

		# ax1 = plt.subplot(122)
		# ax1.set_yticklabels([])
		# #ax.set_xticklabels([])
		# plt2 = plt.contour(x_grid[0,:], y_grid[:,0], measured, levels=clevels, linewidths=0.1)
		# plt.contourf(x_grid[0,:], y_grid[:,0], measured, levels=clevels)
		# plt.clabel(plt2, clevels, colors='k', fontsize=5)
		# plt.title('Measured')
		# #plt.clabel(CS, clevels, inline=True)
		# #plt.contourf(y_grid[:,0], x_grid[0,:], measured, levels=clevels)
		# #plt.axis('equal')
		# #plt.colorbar()
		# plt.xlabel('y/R')
		#plt.ylabel('-x/R')

		#ax = plt.subplot(133)
		cax = fig.add_axes([0.95, 0.11, 0.02, 0.77])
		fig.colorbar(plt1, cax, orientation='vertical', label='SPL [dB]')
		plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name.upper()}_{freq_range.Title}_MR_dBA.pdf', dpi=500, bbox_inches="tight", pad_inches=0.0)
		#plt.show()

def plot_acoustic_contours_tr(plot_name: str, plot_hybrid = False):
	#x_grid, y_grid, measured = read_hart_contour_tecplot(f'{os.path.dirname(os.path.realpath(__file__))}/{plot_name.lower()}-contour-meas.tec')

	namelist = parse_namelist(f'{os.path.dirname(os.path.realpath(__file__))}/{plot_name.upper()}/acoustics/rotor_1/case.nam')
	wopwop_results = parse_wopwop_results_namelist(f'{os.path.dirname(os.path.realpath(__file__))}/{plot_name.upper()}/acoustics/rotor_1', namelist)

	hybrid_wopwop_results = None
	if plot_hybrid:
		hybrid_wopwop_results = parse_wopwop_results_namelist(f'{os.path.dirname(os.path.realpath(__file__))}/../helinovi_hybrid/helinovi/{plot_name.upper()}/acoustics/rotor_1', namelist)

	i_max = len(wopwop_results.oaspl_dba_grid.obs_x)
	j_max = len(wopwop_results.oaspl_dba_grid.obs_x[0])

	phi = np.linspace(0, 2*math.pi, 1000)
	rotor_x = 2*np.cos(phi)
	rotor_z = 2*np.sin(phi)

	for r_idx, freq_range in enumerate(namelist.observers[0].ranges):
		oaspl_linear = [oaspl_db.functions[2].data[0] for oaspl_db in wopwop_results.frequency_ranges_db[r_idx]]

		oaspl_db = [[oaspl_linear[j*i_max + i] for i in range(i_max)] for j in range(j_max)]

		oaspl_linear = [oaspl_db.functions[2].data[0] for oaspl_db in wopwop_results.frequency_ranges_dba[r_idx]]

		oaspl_dba = [[oaspl_linear[j*i_max + i] for i in range(i_max)] for j in range(j_max)]

		y = wopwop_results.oaspl_dba_grid.obs_y[0]
		x = [_x[0] for _x in wopwop_results.oaspl_dba_grid.obs_x]
		x_delta = 4 - x[0]
		x.reverse()
		x = np.asarray(x)
		x = x + x_delta

		clevels = np.linspace(ACOUSTIC_BOUNDS_DB[plot_name][r_idx][0], ACOUSTIC_BOUNDS_DB[plot_name][r_idx][1], ACOUSTIC_BOUNDS_DB[plot_name][r_idx][2])
		cnorm = plt.Normalize(ACOUSTIC_BOUNDS_DB[plot_name][r_idx][0], ACOUSTIC_BOUNDS_DB[plot_name][r_idx][1])

		fig = plt.figure()

		plt2 = plt.contour(y, x, oaspl_db, levels=clevels, linewidths=0.1, cmap=cmap_lines, norm=cnorm)
		plt1 = plt.contourf(y, x, oaspl_db, levels=clevels, cmap=ACOUSTIC_CMAP[plot_name], norm=cnorm)

		plt.ylabel('x/R')
		plt.xlabel('y/R')
		plt.title('Prediction')
		plt.axis('scaled')

		cax = fig.add_axes([0.95, 0.11, 0.02, 0.77])
		fig.colorbar(plt1, cax, orientation='vertical', label='BVI SPL [dB]')
		plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name.upper()}_{freq_range.Title}_TR_dB.pdf', dpi=500, bbox_inches="tight", pad_inches=0.0)

		clevels = np.linspace(ACOUSTIC_BOUNDS_TR_DBA[plot_name][r_idx][0], ACOUSTIC_BOUNDS_TR_DBA[plot_name][r_idx][1], ACOUSTIC_BOUNDS_TR_DBA[plot_name][r_idx][2])
		cnorm = plt.Normalize(ACOUSTIC_BOUNDS_TR_DBA[plot_name][r_idx][0], ACOUSTIC_BOUNDS_TR_DBA[plot_name][r_idx][1])


		num_subplots = 1
		if MEASURED_DATA_TR_DBA[plot_name][r_idx] != "":
			num_subplots = num_subplots + 1

		if UPM_DATA_TR_DBA[plot_name][r_idx] != "" or plot_hybrid:
			num_subplots = num_subplots + 1

		fig = plt.figure()
		_ = plt.subplot(int("1"+str(num_subplots)+"1"))

		plt2 = plt.contour(y, x, oaspl_dba, levels=clevels, linewidths=0.5, cmap=cmap_lines, norm=cnorm)
		plt1 = plt.contourf(y, x, oaspl_dba, levels=clevels, cmap=ACOUSTIC_CMAP[plot_name], norm=cnorm)
		plt.plot(rotor_x, rotor_z, "k", linewidth=1)

		clabels = plt.clabel(plt2, clevels[0::2], colors='k', fontsize=font_size4)
		plt.ylabel('x [m]')
		plt.xlabel('y [m]')
		
		if plot_hybrid:
			plt.title('OpenCOPTER Standard')
		else:
			plt.title('Prediction: OpenCOPTER')

		plt.axis('scaled')
		plt.ylim(-4, 4)
	
		for label in clabels:
			label.set_bbox(dict(facecolor='white', edgecolor='white', pad=0.01))
		
		plotted_subplots = 2

		if MEASURED_DATA_TR_DBA[plot_name][r_idx] != "":
			_ = plt.subplot(int("1"+str(num_subplots)+str(plotted_subplots)))
			arr = plt.imread(f'{os.path.dirname(os.path.realpath(__file__))}/{MEASURED_DATA_TR_DBA[plot_name][r_idx]}')

			plt.imshow(np.flipud(arr) ,interpolation='bilinear', origin='lower', extent=[y[0],y[-1],-4,4])
			plt.plot(rotor_x, rotor_z, "k", linewidth=1)
			plt.title('Measured')
			plt.xlabel('y [m]')
			plotted_subplots = plotted_subplots+1

		if (UPM_DATA_TR_DBA[plot_name][r_idx] != "") and not plot_hybrid:
			_ = plt.subplot(int("1"+str(num_subplots)+str(plotted_subplots)))
			arr = plt.imread(f'{os.path.dirname(os.path.realpath(__file__))}/{UPM_DATA_TR_DBA[plot_name][r_idx]}')
			
			plt.imshow(np.flipud(arr) ,interpolation='bilinear', origin='lower', extent=[y[0],y[-1],-4,4])
			plt.plot(rotor_x, rotor_z, "k", linewidth=1)
			plt.title('Prediction: UPM')
			plt.xlabel('y [m]')
		elif plot_hybrid:
			hybrid_oaspl_linear = [oaspl_db.functions[2].data[0] for oaspl_db in hybrid_wopwop_results.frequency_ranges_dba[r_idx]]

			hybrid_oaspl_dba = [[hybrid_oaspl_linear[j*i_max + i] for i in range(i_max)] for j in range(j_max)]

			_ = plt.subplot(int("1"+str(num_subplots)+str(plotted_subplots)))
			plt2 = plt.contour(y, x, hybrid_oaspl_dba, levels=clevels, linewidths=0.5, cmap=cmap_lines, norm=cnorm)
			plt1 = plt.contourf(y, x, hybrid_oaspl_dba, levels=clevels, cmap=ACOUSTIC_CMAP[plot_name], norm=cnorm)
			plt.plot(rotor_x, rotor_z, "k", linewidth=1)

			clabels = plt.clabel(plt2, clevels[0::2], colors='k', fontsize=font_size4)
			for label in clabels:
				label.set_bbox(dict(facecolor='white', edgecolor='white', pad=0.01))

			plt.xlabel('y [m]')
			
			plt.title('OpenCOPTER Hybrid')
			plt.axis('scaled')
			plt.ylim(-4, 4)

		cax = fig.add_axes([0.95, 0.31, 0.02, 0.37])
		fig.colorbar(plt1, cax, orientation='vertical', label='SPL [dB]')
		if not plot_hybrid:
			plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name.upper()}_{freq_range.Title}_TR_dBA.pdf', dpi=500, bbox_inches="tight", pad_inches=0.0)
		else:
			plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name.upper()}_{freq_range.Title}_TR_hybrid_dBA.pdf', dpi=500, bbox_inches="tight", pad_inches=0.0)

		plt.clf()
		plt.cla()

SPAN_LOCATION = [
	0.80,
	0.87,
	0.97
]

ROTOR = [
	"MR",
	"TR"
]

MEASURED_LOADING = {
	"ID1": [
		["xx", "ID1_MR_loading_87.tsv", "xx"],
		["ID1_TR_loading_80_5r.tsv", "xx", "xx"]
	],
	"ID1_ASD": [
		["xx", "xx", "xx"],
		["xx", "xx", "ID1_ASD_TR_loading_97_5r.tsv"]
	],
	"ID1_ASU": [
		["xx", "xx", "xx"],
		["xx", "xx", "ID1_ASU_TR_loading_97_5r.tsv"]
	],
	"ID2": [
		["xx", "ID2_MR_loading_87.tsv", "xx"],
		["ID2_TR_loading_80_5r.tsv", "xx", "ID2_TR_loading_97_5r.tsv"]
	],
	"ID2_ASD": [
		["xx", "xx", "xx"],
		["xx", "xx", "xx"]
	],
	"ID2_ASU": [
		["xx", "xx", "xx"],
		["xx", "xx", "xx"]
	],
	"ID5": [
		["xx", "ID5_MR_loading_87.tsv", "xx"],
		["xx", "xx", "xx"]
	]
}

MEASURED_LOADING_UPM = {
	"ID1": [
		["xx", "xx", "xx"],
		["ID1_TR_loading_80_5r_UPM.tsv", "xx", "xx"]
	],
	"ID1_ASD": [
		["xx", "xx", "xx"],
		["xx", "xx", "ID1_ASD_TR_loading_97_5r_UPM.tsv"]
	],
	"ID1_ASU": [
		["xx", "xx", "xx"],
		["xx", "xx", "ID1_ASU_TR_loading_97_5r_UPM.tsv"]
	],
	"ID2": [
		["xx", "xx", "xx"],
		["ID2_TR_loading_80_5r_UPM.tsv", "xx", "ID2_TR_loading_97_5r_UPM.tsv"]
	],
	"ID2_ASD": [
		["xx", "xx", "xx"],
		["xx", "xx", "xx"]
	],
	"ID2_ASU": [
		["xx", "xx", "xx"],
		["xx", "xx", "xx"]
	],
	"ID5": [
		["xx", "xx", "xx"],
		["xx", "xx", "xx"]
	]
}

MEASURED_LOADING_TR = {
	"ID1": ["ID1_TR_loading_80_1r.tsv", "xx", "ID1_TR_loading_97_1r.tsv"],
	"ID1_ASU": ["xx", "xx", "ID1_ASU_TR_loading_97_1r.tsv"],
	"ID1_ASD": ["xx", "xx", "ID1_ASD_TR_loading_97_1r.tsv"],
	"ID2": ["ID2_TR_loading_80_1r.tsv", "xx", "ID2_TR_loading_97_1r.tsv"],
	"ID2_ASU": ["xx", "xx", "xx"],
	"ID2_ASD": ["xx", "xx", "xx"],
	"ID5": ["xx", "xx", "xx"]
}

MEASURED_LOADING_TR_UPM = {
	"ID1": ["ID1_TR_loading_80_1r_UPM.tsv", "xx", "ID1_TR_loading_97_1r_UPM.tsv"],
	"ID1_ASU": ["xx", "xx", "ID1_ASU_TR_loading_97_1r_UPM.tsv"],
	"ID1_ASD": ["xx", "xx", "ID1_ASD_TR_loading_97_1r_UPM.tsv"],
	"ID2": ["ID2_TR_loading_80_1r_UPM.tsv", "xx", "ID2_TR_loading_97_1r_UPM.tsv"],
	"ID2_ASU": ["xx", "xx", "xx"],
	"ID2_ASD": ["xx", "xx", "xx"],
	"ID5": ["xx", "xx", "xx"]
}

MEASURED_WAKE = {
	"ID1": ["HeliNOVI_ID1_MR_wake_trajectory_plane0_measured.tsv", "HeliNOVI_ID1_MR_wake_trajectory_plane1_measured.tsv"],
	"ID1_ASU": ["xx", "xx"],
	"ID1_ASD": ["xx", "xx"],
	"ID2": ["HeliNOVI_ID2_MR_wake_trajectory_plane0_measured.tsv", "HeliNOVI_ID2_MR_wake_trajectory_plane1_measured.tsv"],
	"ID2_ASU": ["xx", "xx"],
	"ID2_ASD": ["xx", "xx"],
	"ID5": ["xx", "xx"]
}

def plot_mr_wake(plot_name: str):
	#with open(f'{os.path.dirname(os.path.realpath(__file__))}/../helinovi_params.json5') as param_file:
	#	params = pyjson5.load(param_file)

	#flight_condition = list(filter(lambda x: x["name"] == plot_name, params['flight_conditions']))[0]

	blade_results = loadmat(f'{os.path.dirname(os.path.realpath(__file__))}/{plot_name.upper()}/results.mat')

	wake_element_piv = blade_results['wake_element_piv']

	# num_slices = wake_x.shape[0]
	# num_samples = wake_x.shape[1]

	phi = np.linspace(0, 2*math.pi, 1000)
	rotor_x = 0.383*np.cos(phi)
	rotor_z = 0.383*np.sin(phi)

	# print(f"num_slices: {num_slices}")
	# print(f"num_samples: {num_samples}")
	markers = ['x', '+']
	line_colors = ['b', 'r']

	plt.figure(num = 1)
	
	for b_idx in range(4):

		plt.plot(rotor_x, rotor_z, 'k', label='_nolegend_', linewidth=1.0)

		for s_idx in range(2):
		
			wake_x = -wake_element_piv[0, b_idx, s_idx, 0, :]
			wake_z = wake_element_piv[0, b_idx, s_idx, 1, :]

			in_wake_x = np.logical_and(wake_x <= 0.4, wake_x >= -0.4)
			in_wake_z = np.logical_and(wake_z <= 0.4, wake_z >= -0.4)

			non_zero_x = wake_x != 0.0
			non_zero_z = wake_z != 0.0
			non_zero = np.logical_and(non_zero_z, non_zero_x)

			in_box = np.logical_and(in_wake_x, in_wake_z)

			good_wake_points = np.logical_and(non_zero, in_box)
			
			wake_x = wake_x[good_wake_points]
			wake_z = wake_z[good_wake_points]

			if os.path.exists(f'{os.path.dirname(os.path.realpath(__file__))}/{MEASURED_WAKE[plot_name][s_idx]}'):
				measured_wake_data = read_hart_blade_data_tsv(f'{os.path.dirname(os.path.realpath(__file__))}/{MEASURED_WAKE[plot_name][s_idx]}')
				measured_x = measured_wake_data[0,:]
				measured_z = measured_wake_data[1,:]

			plt.plot(wake_x[:], wake_z[:], line_colors[s_idx], linewidth=0.5, markersize=5.0, label=f"OpenCOPTER Plane {s_idx}")
			plt.plot(measured_x, measured_z, f"{line_colors[s_idx]}{markers[s_idx]}-", linewidth=0.5, markersize=5.0, label=f"Measured Plane {s_idx}")

		plt.axis('square')	
		plt.xlim(-.5, .5)
		plt.ylim(-.5, .5)
		
		plt.xlabel("x [m]")
		plt.ylabel("z [m]")
		plt.legend(ncol=2, frameon=False, fancybox=False, fontsize=font_size4+2)
		#plt.legend(["OpenCOPTER", "Measured"])

		plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name}_{ROTOR[0]}_wake_trajectory_blade_{b_idx}.pdf', dpi=500, bbox_inches="tight", pad_inches=0.0)
		plt.cla()
		plt.clf()

def plot_hybrid_wake(plot_name: str, azimuth: int = 0):
	blade_results = loadmat(f'{os.path.dirname(os.path.realpath(__file__))}/{plot_name.upper()}/results.mat')
	hybrid_blade_results = loadmat(f'{os.path.dirname(os.path.realpath(__file__))}/../helinovi_hybrid/helinovi/{plot_name.upper()}/results.mat')

	b1_tip_x = -hybrid_blade_results["rotor_1_wake_timehistory"][azimuth, 0, 0, 0]
	b1_tip_y = hybrid_blade_results["rotor_1_wake_timehistory"][azimuth, 0, 1, 0]
	
	azi1 = math.atan2(b1_tip_y, b1_tip_x)
	
	b1_root_y = 0.383*0.42*math.sin(azi1)
	b1_root_x = 0.383*0.42*math.cos(azi1)
	
	b2_tip_x = -hybrid_blade_results["rotor_1_wake_timehistory"][azimuth, 1, 0, 0]
	b2_tip_y = hybrid_blade_results["rotor_1_wake_timehistory"][azimuth, 1, 1, 0]
	
	azi2 = math.atan2(b2_tip_y, b2_tip_x)
	b2_root_y = 0.383*0.41*math.sin(azi2)
	b2_root_x = 0.383*0.41*math.cos(azi2)
	
	plt.plot([b2_root_x, b2_tip_x], [b2_root_y, b2_tip_y], 'k', linewidth=2.0)

	plt.plot([b1_root_x, b1_tip_x], [b1_root_y, b1_tip_y], 'k', linewidth=2.0)
	plt.plot(-hybrid_blade_results['rotor_1_wake_timehistory'][azimuth, 0, 0, :], hybrid_blade_results["rotor_1_wake_timehistory"][azimuth, 0, 1, :], 'r', linewidth=1.0, label="Hybrid blade 0")
	plt.plot(-hybrid_blade_results['rotor_1_wake_timehistory'][azimuth, 1, 0, :], hybrid_blade_results["rotor_1_wake_timehistory"][azimuth, 1, 1, :], 'b', linewidth=1.0, label="Hybrid blade 1")
	
	plt.plot(-blade_results["rotor_1_wake_timehistory"][azimuth, 0, 0, :], blade_results["rotor_1_wake_timehistory"][azimuth, 0, 1, :], 'r--', linewidth=1.0, label="Standard blade 0")
	plt.plot(-blade_results["rotor_1_wake_timehistory"][azimuth, 1, 0, :], blade_results["rotor_1_wake_timehistory"][azimuth, 1, 1, :], 'b--', linewidth=1.0, label="Standard blade 1")
	
	#plt.title(f"iteration %d\n", azimuth))
	if plot_name == "ID1":
		plt.legend(ncol=2, frameon=False, fancybox=False, fontsize=font_size4 + 4)

	plt.axis('square')
	plt.xlim([-0.5, 0.5])
	plt.ylim([-0.5, 0.5])
	plt.ylabel('$z$ [m]')
	plt.ylabel('$x$ [m]')
	plt.xticks(fontsize=font_size3)
	plt.yticks(fontsize=font_size3)
	plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name}_TR_hybrid_wake_compare.pdf', dpi=500, bbox_inches="tight", pad_inches=0.0)
	plt.cla()
	plt.clf()

def plot_blade_normal_pressures(plot_name: str, hybrid_compare = False, isolated_compare = False):

	with open(f'{os.path.dirname(os.path.realpath(__file__))}/../helinovi_params.json5') as param_file:
		params = pyjson5.load(param_file)

	flight_condition = list(filter(lambda x: x["name"] == plot_name, params['flight_conditions']))[0]

	omegas = []
	for m in flight_condition["motion"]:
		if ('blade_element_func' in m) and (m['blade_element_func'] == "rotation"):
			omegas.append(m['omega'])
	
	blade_results = loadmat(f'{os.path.dirname(os.path.realpath(__file__))}/{plot_name.upper()}/results.mat')

	hybrid_blade_results = None
	if hybrid_compare:
		hybrid_blade_results = loadmat(f'{os.path.dirname(os.path.realpath(__file__))}/../helinovi_hybrid/helinovi/{plot_name.upper()}/results.mat')

	isolated_blade_results = None
	if isolated_compare:
		isolated_blade_results = loadmat(f'{os.path.dirname(os.path.realpath(__file__))}/../helinovi_tr/{plot_name.upper()}/results.mat')
	
	plt.figure(num=1)

	for rotor_idx, rotor_blade_loading in enumerate(blade_results['span_element_af_loading']):
		for span_idx, blade_loading in enumerate(rotor_blade_loading):

			omega = omegas[rotor_idx]

			measured_blade_loading = np.asarray([])#np.zeros(1440)
			measured_azimuth = np.asarray([])#np.linspace(0, 360, measured_blade_loading.size)

			plots_done = 1

			computational_azimuth = blade_results['span_element_azimuth'][rotor_idx]

			delta_azi = 0

			if os.path.exists(f'{os.path.dirname(os.path.realpath(__file__))}/{MEASURED_LOADING_UPM[plot_name][rotor_idx][span_idx]}'):
				upm_loading = read_hart_blade_data_tsv(f'{os.path.dirname(os.path.realpath(__file__))}/{MEASURED_LOADING_UPM[plot_name][rotor_idx][span_idx]}')
				upm_azimuth = upm_loading[0,:]
				upm_blade_loading = upm_loading[1,:]

				plt.plot(upm_azimuth, upm_blade_loading - upm_blade_loading.mean(), 'ks-', linewidth=0.5, markersize=0.5, mfc='none', label='UPM')
				plots_done = plots_done + 1

			if os.path.exists(f'{os.path.dirname(os.path.realpath(__file__))}/{MEASURED_LOADING[plot_name][rotor_idx][span_idx]}'):
				measured_loading = read_hart_blade_data_tsv(f'{os.path.dirname(os.path.realpath(__file__))}/{MEASURED_LOADING[plot_name][rotor_idx][span_idx]}')
				measured_azimuth = measured_loading[0,:]
				measured_blade_loading = measured_loading[1,:]

				plots_done = plots_done + 1

				plt.plot(measured_azimuth, measured_blade_loading - measured_blade_loading.mean(), 'ro-', linewidth=0.5, markersize=1.0, mfc='none', label='Measured')

			plt.plot(computational_azimuth + delta_azi, blade_loading - blade_loading.mean(), 'b-', linewidth=0.5, markersize=1.0, mfc='none', label="OpenCOPTER")
			
			if rotor_idx == 0:
				plt.legend(ncol=plots_done, frameon=False, fancybox=False, fontsize=font_size4 + 4)
				plt.ylim(-0.6, 0.7)
			else:
				plt.legend(ncol=plots_done, frameon=False, fancybox=False, fontsize=font_size4 + 2)
				plt.ylim(-0.5, 0.5)

			plt.xlabel("$\psi$ [$^\circ$]")
			plt.ylabel("$C_N$")
			if rotor_idx == 1:
				plt.gca().set_aspect(700)
				plt.xticks(ticks=np.linspace(0, 1800, 6))
			else:
				plt.gca().set_aspect(200)
				plt.xticks(ticks=np.linspace(0, 360, 5))

			plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name}_{ROTOR[rotor_idx]}_{SPAN_LOCATION[span_idx]:0.2f}_normal_loading.pdf', dpi=500, bbox_inches="tight", pad_inches=0.0)
			plt.cla()
			plt.clf()

			if hybrid_compare and rotor_idx == 1:
				hybrid_loading = hybrid_blade_results['span_element_af_loading'][rotor_idx][span_idx]
				hybrid_azimuth = hybrid_blade_results['span_element_azimuth'][rotor_idx]

				#plt.figure(num=1)
				plots_done = 2
				if os.path.exists(f'{os.path.dirname(os.path.realpath(__file__))}/{MEASURED_LOADING[plot_name][rotor_idx][span_idx]}'):
					plots_done = plots_done + 1

					plt.plot(measured_azimuth, measured_blade_loading - measured_blade_loading.mean(), 'ro-', linewidth=0.5, markersize=1.0, mfc='none', label='Measured')

				plt.plot(
					computational_azimuth + delta_azi, blade_loading - blade_loading.mean(), 'b-',
					linewidth=0.5, markersize=1.0, mfc='none',
					label='OpenCOPTER Standard'
				)
				plt.plot(
					hybrid_azimuth + delta_azi, hybrid_loading - hybrid_loading.mean(), 'k-',
					linewidth=0.5, markersize=1.0, mfc='none',
					label='OpenCOPTER Hybrid'
					
				)

				#plt.title(f'{TITLE_DICT[plot_name.upper()]} {SPAN_LOCATION[span_idx]*100}% Span blade loading. Mean removed.')
				plt.legend(ncol=plots_done, frameon=False, fancybox=False, fontsize=font_size4 + 2)
				plt.ylim(-0.5, 0.5)
				plt.gca().set_aspect(700)
				plt.xlabel("$\psi$ [$^\circ$]")
				plt.ylabel("$C_N$")
				plt.xticks(ticks=np.linspace(0, 1800, 6))
				plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name}_{ROTOR[rotor_idx]}_{SPAN_LOCATION[span_idx]:0.2f}_normal_loading_hybrid_compare.pdf', dpi=500, bbox_inches="tight", pad_inches=0.0)
				plt.cla()
				plt.clf()

			if isolated_compare and rotor_idx == 1:
				isolated_loading = isolated_blade_results['span_element_af_loading'][span_idx]
				isolated_azimuth = isolated_blade_results['span_element_azimuth']

				#plt.figure(num=1)
				plots_done = 2
				# if os.path.exists(f'{os.path.dirname(os.path.realpath(__file__))}/{MEASURED_LOADING[plot_name][rotor_idx][span_idx]}'):
				# 	plots_done = plots_done + 1

				# 	plt.plot(measured_azimuth, measured_blade_loading - measured_blade_loading.mean(), 'ro-', linewidth=0.5, markersize=1.0, mfc='none', label='Measured')

				plt.plot(
					computational_azimuth + delta_azi, blade_loading - blade_loading.mean(), 'b-',
					linewidth=0.5, markersize=1.0, mfc='none',
					label='OpenCOPTER MR+TR'
				)
				plt.plot(
					isolated_azimuth + delta_azi, isolated_loading - isolated_loading.mean(), 'r-',
					linewidth=0.5, markersize=1.0, mfc='none',
					label='OpenCOPTER Isolated TR'
					
				)

				#plt.title(f'{TITLE_DICT[plot_name.upper()]} {SPAN_LOCATION[span_idx]*100}% Span blade loading. Mean removed.')
				plt.legend(ncol=plots_done, frameon=False, fancybox=False, fontsize=font_size4 + 2)
				plt.ylim(-0.5, 0.5)
				plt.gca().set_aspect(700)
				plt.xlabel("$\psi$ [$^\circ$]")
				plt.ylabel("$C_N$")
				plt.xticks(ticks=np.linspace(0, 1800, 6))
				plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name}_{ROTOR[rotor_idx]}_{SPAN_LOCATION[span_idx]:0.2f}_normal_loading_isolated_compare.pdf', dpi=500, bbox_inches="tight", pad_inches=0.0)
				plt.cla()
				plt.clf()

			if rotor_idx == 1:
				start_azi = 360 - delta_azi
				measured_azimuth_tr = []
				measured_blade_loading_tr = []

				#plt.figure(num=1)

				plots_done = 1
				if os.path.exists(f'{os.path.dirname(os.path.realpath(__file__))}/{MEASURED_LOADING_TR_UPM[plot_name][span_idx]}'):
					upm_loading_tr = read_hart_blade_data_tsv(f'{os.path.dirname(os.path.realpath(__file__))}/{MEASURED_LOADING_TR_UPM[plot_name][span_idx]}')
					upm_azimuth_tr = upm_loading_tr[0,:]
					upm_blade_loading_tr = upm_loading_tr[1,:]
				
					plt.plot(upm_azimuth_tr, upm_blade_loading_tr - np.asarray(upm_blade_loading_tr).mean(), "ks-", linewidth=0.5, markersize=0.5, label='UPM')
					plots_done = plots_done + 1

				if os.path.exists(f'{os.path.dirname(os.path.realpath(__file__))}/{MEASURED_LOADING_TR[plot_name][span_idx]}'):
					measured_loading_tr = read_hart_blade_data_tsv(f'{os.path.dirname(os.path.realpath(__file__))}/{MEASURED_LOADING_TR[plot_name][span_idx]}')
					measured_azimuth_tr = measured_loading_tr[0,:]
					measured_blade_loading_tr = measured_loading_tr[1,:]
				
					plt.plot(measured_azimuth_tr, measured_blade_loading_tr - np.asarray(measured_blade_loading_tr).mean(), "r.-", linewidth=0.5, markersize=1.5, label='Measured')
					plots_done = plots_done + 1

				plt.plot(computational_azimuth[0:360], blade_loading[start_azi:start_azi+360] - np.asarray(blade_loading[start_azi:start_azi+360]).mean(), 'b', linewidth=0.5, label='OpenCOPTER')
				#plt.title(f'{TITLE_DICT[plot_name.upper()]} {SPAN_LOCATION[span_idx]*100}% Span blade loading. Mean removed.')
				plt.legend(ncol=plots_done, frameon=False, fancybox=False, fontsize=font_size4 + 2)

				plt.ylim(-0.5, 0.5)
				plt.xticks(ticks=np.linspace(0, 360, 5))
				plt.gca().set_aspect(120)
				plt.xlabel("$\psi$ [$^\circ$]")
				plt.ylabel("$C_N$")
				plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name}_{ROTOR[rotor_idx]}_{SPAN_LOCATION[span_idx]:0.2f}_normal_loading_1rev.pdf', dpi=500, bbox_inches="tight", pad_inches=0.0)
				plt.cla()
				plt.clf()

				if hybrid_compare:
					hybrid_loading = hybrid_blade_results['span_element_af_loading'][rotor_idx][span_idx]
					hybrid_azimuth = hybrid_blade_results['span_element_azimuth'][rotor_idx]

					#plt.figure(num=1)
					plots_done = 2
					if os.path.exists(f'{os.path.dirname(os.path.realpath(__file__))}/{MEASURED_LOADING_TR[plot_name][span_idx]}'):
						plt.plot(measured_azimuth_tr, measured_blade_loading_tr - np.asarray(measured_blade_loading_tr).mean(), "r.-", linewidth=0.5, markersize=1.0, label='Measured')
						plots_done = plots_done + 1

					plt.plot(
						computational_azimuth[0:360], blade_loading[start_azi:start_azi+360] - np.asarray(blade_loading[start_azi:start_azi+360]).mean(), 'b',
						linewidth=0.5,
						label='OpenCOPTER Standard'
					)
					plt.plot(
						hybrid_azimuth[0:360], hybrid_loading[start_azi:start_azi+360] - np.asarray(hybrid_loading[start_azi:start_azi+360]).mean(), 'k',
						linewidth=0.5,
						label='OpenCOPTER Hybrid'
					)
					
					#plt.title(f'{TITLE_DICT[plot_name.upper()]} {SPAN_LOCATION[span_idx]*100}% Span blade loading. Mean removed.')
					plt.legend(ncol=plots_done, frameon=False, fancybox=False, fontsize=font_size4 + 2)
					plt.ylim(-0.5, 0.5)
					plt.xticks(ticks=np.linspace(0, 360, 5))
					plt.gca().set_aspect(120)
					plt.xlabel("$\psi$ [$^\circ$]")
					plt.ylabel("$C_N$")
					plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name}_{ROTOR[rotor_idx]}_{SPAN_LOCATION[span_idx]:0.2f}_normal_loading_1rev_hybrid_compare.pdf', dpi=500, bbox_inches="tight", pad_inches=0.0)
					plt.cla()
					plt.clf()

				if isolated_compare:
					isolated_loading = isolated_blade_results['span_element_af_loading'][span_idx]
					isolated_azimuth = isolated_blade_results['span_element_azimuth']

					#plt.figure(num=1)
					plots_done = 2
					# if os.path.exists(f'{os.path.dirname(os.path.realpath(__file__))}/{MEASURED_LOADING_TR[plot_name][span_idx]}'):
					# 	plt.plot(measured_azimuth_tr, measured_blade_loading_tr - np.asarray(measured_blade_loading_tr).mean(), "r.-", linewidth=0.5, label='Measured')
					# 	plots_done = plots_done + 1

					plt.plot(
						computational_azimuth[0:360], blade_loading[start_azi:start_azi+360] - np.asarray(blade_loading[start_azi:start_azi+360]).mean(), 'b',
						linewidth=0.5,
						label='OpenCOPTER MR+TR'
					)
					plt.plot(
						isolated_azimuth[0:360], isolated_loading[start_azi:start_azi+360] - np.asarray(isolated_loading[start_azi:start_azi+360]).mean(), 'r',
						linewidth=0.5,
						label='OpenCOPTER Isolated TR'
					)
					
					#plt.title(f'{TITLE_DICT[plot_name.upper()]} {SPAN_LOCATION[span_idx]*100}% Span blade loading. Mean removed.')
					plt.legend(ncol=plots_done, frameon=False, fancybox=False, fontsize=font_size4 + 2)
					plt.ylim(-0.5, 0.5)
					plt.xticks(ticks=np.linspace(0, 360, 5))
					plt.gca().set_aspect(120)
					plt.xlabel("$\psi$ [$^\circ$]")
					plt.ylabel("$C_N$")
					plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name}_{ROTOR[rotor_idx]}_{SPAN_LOCATION[span_idx]:0.2f}_normal_loading_1rev_isolated_compare.pdf', dpi=500, bbox_inches="tight", pad_inches=0.0)
					plt.cla()
					plt.clf()
				
			# 10/rev high pass
			if measured_blade_loading.shape[0] != 0:
				measured_dt = (2.0*math.pi/omega)/measured_blade_loading.shape[0]
				computed_dt = (2.0*math.pi/omega)/blade_loading.shape[0]

				print(f'measured_dt: {measured_dt}')
				print(f'computed_dt: {computed_dt}')

				fblade_loading = sft.fft(blade_loading - blade_loading.mean())
				fblade_freqz = sft.fftfreq(fblade_loading.shape[0], d=computed_dt)

				fmeasured_blade_loading = sft.fft(measured_blade_loading - measured_blade_loading.mean())
				fmeasured_blade_freqz = sft.fftfreq(fmeasured_blade_loading.shape[0], d=measured_dt)

				ten_per_rev = 10*omega*RADPS_2_HZ

				#print(f'fblade_freqz: {fblade_freqz}')
				#print(f'fmeasured_blade_freqz: {fmeasured_blade_freqz}')

				measured_sos = sig.butter(13, ten_per_rev, 'highpass', output='sos', fs=1/measured_dt)
				computed_sos = sig.butter(13, ten_per_rev, 'highpass', output='sos', fs=1/computed_dt)

				blade_loading = sig.sosfilt(computed_sos, blade_loading - blade_loading.mean())
				measured_blade_loading = sig.sosfilt(measured_sos, measured_blade_loading - measured_blade_loading.mean())

			#plt.figure(num=1)
			plt.plot(measured_azimuth, measured_blade_loading, computational_azimuth, blade_loading, linewidth=0.5)
			plt.title(f'{TITLE_DICT[plot_name.upper()]}: 87% Span blade loading. 10/rev highpassed.')
			plt.legend(['Measured', 'OpenCOPTER'])
			#plt.xlim([20, 90])
			plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name}_{ROTOR[rotor_idx]}_{SPAN_LOCATION[span_idx]:0.2f}_normal_loading_hp.pdf', dpi=500, bbox_inches="tight", pad_inches=0.0)
			plt.cla()
			plt.clf()

			aoa_eff = blade_results['span_element_aoa_eff'][rotor_idx][span_idx]
			
			#plt.figure(num=1)
			plt.plot(computational_azimuth, aoa_eff, linewidth=0.5)
			plt.title(f'{TITLE_DICT[plot_name.upper()]} 87% Effective AOA.')
			#plt.legend(['Measured', 'OpenCOPTER'])
			plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name}_{ROTOR[rotor_idx]}_{SPAN_LOCATION[span_idx]:0.2f}_aoa_eff.pdf', dpi=500, bbox_inches="tight", pad_inches=0.0)
			plt.cla()
			plt.clf()

			u_p = blade_results['span_element_up'][rotor_idx][span_idx]
			
			#plt.figure(num=1)
			plt.plot(computational_azimuth, u_p, linewidth=0.5)
			plt.title(f'{TITLE_DICT[plot_name.upper()]} 87% U_p.')
			#plt.legend(['Measured', 'OpenCOPTER'])
			plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name}_{ROTOR[rotor_idx]}_{SPAN_LOCATION[span_idx]:0.2f}_u_p.pdf', dpi=500, bbox_inches="tight", pad_inches=0.0)
			plt.cla()
			plt.clf()

			aoa = blade_results['span_element_aoa'][rotor_idx][span_idx]
			
			#plt.figure(num=1)
			plt.plot(computational_azimuth, aoa, linewidth=0.5)
			plt.title(f'{TITLE_DICT[plot_name.upper()]} 87% AOA.')
			#plt.legend(['Measured', 'OpenCOPTER'])
			plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name}_{ROTOR[rotor_idx]}_{SPAN_LOCATION[span_idx]:0.2f}_aoa.pdf', dpi=500, bbox_inches="tight", pad_inches=0.0)
			plt.cla()
			plt.clf()

			theta = blade_results['span_element_theta'][rotor_idx][span_idx]
			
			#plt.figure(num=1)
			plt.plot(computational_azimuth, theta, linewidth=0.5)
			plt.title(f'{TITLE_DICT[plot_name.upper()]} 87% theta.')
			#plt.legend(['Measured', 'OpenCOPTER'])
			plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name}_{ROTOR[rotor_idx]}_{SPAN_LOCATION[span_idx]:0.2f}_theta.pdf', dpi=500, bbox_inches="tight", pad_inches=0.0)
			plt.cla()
			plt.clf()

			inflow_angle = blade_results['span_element_inflow_angle'][rotor_idx][span_idx]
			
			#plt.figure(num=1)
			plt.plot(computational_azimuth, inflow_angle, linewidth=0.5)
			plt.title(f'{TITLE_DICT[plot_name.upper()]} 87% inflow angle.')
			#plt.legend(['Measured', 'OpenCOPTER'])
			plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name}_{ROTOR[rotor_idx]}_{SPAN_LOCATION[span_idx]:0.2f}_inflow_angle.pdf', dpi=500, bbox_inches="tight", pad_inches=0.0)
			plt.cla()
			plt.clf()

			gamma = blade_results['span_element_gamma'][rotor_idx][span_idx]
			
			#plt.figure(num=1)
			plt.plot(computational_azimuth, gamma, linewidth=0.5)
			plt.title(f'{TITLE_DICT[plot_name.upper()]} 87% gamma.')
			#plt.legend(['Measured', 'OpenCOPTER'])
			plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name}_{ROTOR[rotor_idx]}_{SPAN_LOCATION[span_idx]:0.2f}_gamma.pdf', dpi=500, bbox_inches="tight", pad_inches=0.0)
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

		#plt.figure(num=1)
		plt.plot(measured_x, measured_z, '*', x_tppc, z_tppc, '-', linewidth=0.5)
		plt.title(f'{TITLE_DICT[plot_name.upper()]}: Y slice: {target_y_slice:.2f}')
		plt.legend(['Measured', 'OpenCOPTER'])
		#plt.axis('scaled')
		plt.xlim([-1.1, 1.1])
		plt.ylim([-0.1, 0.15])
		plt.gca().invert_xaxis()
		plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/Hart_{plot_name}_wake_y{target_y_slice:.2f}.pdf', dpi=500, bbox_inches="tight", pad_inches=0.0)
		plt.cla()
		plt.clf()

		core_size = wake_element_core_size[idx, 0:wake_element_index[idx]]
		#plt.figure(num=1)
		plt.plot(core_size, linewidth=0.5)
		plt.title(f'{TITLE_DICT[plot_name.upper()]}: Core size, Y slice: {target_y_slice:.2f}')
		#plt.legend(['Measured', 'OpenCOPTER'])
		#plt.axis('scaled')
		#plt.xlim([-1.1, 1.1])
		#plt.ylim([-0.1, 0.15])
		#plt.gca().invert_xaxis()
		plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/Hart_{plot_name}_core_size_y{target_y_slice:.2f}.pdf', dpi=500, bbox_inches="tight", pad_inches=0.0)
		plt.cla()
		plt.clf()

	wake_fil_core_size = wake_results['wake_0_core_size'][2,:]
	#core_size = wake_element_core_size[idx, 0:wake_element_index[idx]]
	#plt.figure(num=1)
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

	with open(f'{os.path.dirname(os.path.realpath(__file__))}/../helinovi_params.json5') as param_file:
		params = pyjson5.load(param_file)

	flight_condition = list(filter(lambda x: x["name"] == plot_name, params['flight_conditions']))[0]

	omegas = []
	for m in flight_condition["motion"]:
		if ('blade_element_func' in m) and (m['blade_element_func'] == "rotation"):
			omegas.append(m['omega'])

	omega = omegas[0]
	ten_per_rev = 10*omega*RADPS_2_HZ

	blade_results = loadmat(f'{os.path.dirname(os.path.realpath(__file__))}/{plot_name.upper()}/results.mat')

	for rotor_idx, _ in enumerate(blade_results['blade_twist_array']):
		blade_idx = 2
		if rotor_idx == 1:
			blade_idx = 0

		blade_twist_array = blade_results['blade_twist_array'][rotor_idx, blade_idx,:]
		blade_twist_azimuth = blade_results['blade_twist_azimuth'][rotor_idx, blade_idx,:]
		

		computed_dt = (2.0*math.pi/omega)/blade_twist_array.shape[0]
		computed_sos = sig.butter(13, ten_per_rev, 'highpass', output='sos', fs=1/computed_dt)

		#print(f'blade_twist_array: {blade_twist_array}')
		plt.figure(num=1)
		plt.plot(blade_twist_azimuth, blade_twist_array*(180/math.pi), linewidth=0.5)
		plt.title(f'{TITLE_DICT[plot_name.upper()]}: Blade twist over time.')
		plt.xlim([0, 360])
		plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name}_{ROTOR[rotor_idx]}_blade_twist.pdf', dpi=500, bbox_inches="tight", pad_inches=0.0)
		plt.cla()
		plt.clf()
		
		elastic_twist_array = blade_results['elastic_twist_array'][rotor_idx, blade_idx,:]
		
		#print(f'blade_twist_array: {blade_twist_array}')
		plt.figure(num=1)
		plt.plot(blade_twist_azimuth, elastic_twist_array, linewidth=0.5)
		plt.title(f'{TITLE_DICT[plot_name.upper()]}: Elastic twist over time.')
		plt.xlim([0, 360])
		plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name}_{ROTOR[rotor_idx]}_elastic_twist.pdf', dpi=500, bbox_inches="tight", pad_inches=0.0)
		plt.cla()
		plt.clf()

		collective_pitch_array = blade_results['collective_pitch_array'][rotor_idx, :]
		
		#print(f'blade_twist_array: {blade_twist_array}')
		plt.figure(num=1)
		plt.plot(blade_twist_azimuth, collective_pitch_array, linewidth=0.5)
		plt.title(f'{TITLE_DICT[plot_name.upper()]}: Collective pitch over time.')
		plt.xlim([0, 360])
		plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name}_{ROTOR[rotor_idx]}_collective.pdf', dpi=500, bbox_inches="tight", pad_inches=0.0)
		plt.cla()
		plt.clf()

		sin_pitch_array = blade_results['sin_pitch_array'][rotor_idx, blade_idx,:]
		sin_pitch_array_hp = sig.sosfilt(computed_sos, sin_pitch_array)

		#print(f'blade_twist_array: {blade_twist_array}')
		plt.figure(num=1)
		plt.plot(blade_twist_azimuth, sin_pitch_array, linewidth=0.5)
		plt.title(f'{TITLE_DICT[plot_name.upper()]}: 1S over time.')
		plt.xlim([0, 360])
		plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name}_{ROTOR[rotor_idx]}_1s.pdf', dpi=500, bbox_inches="tight", pad_inches=0.0)
		plt.cla()
		plt.clf()

		#print(f'blade_twist_array: {blade_twist_array}')
		plt.figure(num=1)
		plt.plot(blade_twist_azimuth, sin_pitch_array_hp, linewidth=0.5)
		plt.title(f'{TITLE_DICT[plot_name.upper()]}: 1S over time HP.')
		plt.xlim([0, 360])
		plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name}_{ROTOR[rotor_idx]}_1s_hp.pdf', dpi=500, bbox_inches="tight", pad_inches=0.0)
		plt.cla()
		plt.clf()

		cos_pitch_array = blade_results['cos_pitch_array'][rotor_idx, blade_idx,:]
		cos_pitch_array_hp = sig.sosfilt(computed_sos, cos_pitch_array)

		#print(f'blade_twist_array: {blade_twist_array}')
		plt.figure(num=1)
		plt.plot(blade_twist_azimuth, cos_pitch_array, linewidth=0.5)
		plt.title(f'{TITLE_DICT[plot_name.upper()]}: 1C over time.')
		plt.xlim([0, 360])
		plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name}_{ROTOR[rotor_idx]}_1c.pdf', dpi=500, bbox_inches="tight", pad_inches=0.0)
		plt.cla()
		plt.clf()

		plt.figure(num=1)
		plt.plot(blade_twist_azimuth, cos_pitch_array_hp, linewidth=0.5)
		plt.title(f'{TITLE_DICT[plot_name.upper()]}: 1C over time HP.')
		plt.xlim([0, 360])
		plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name}_{ROTOR[rotor_idx]}_1c_hp.pdf', dpi=500, bbox_inches="tight", pad_inches=0.0)
		plt.cla()
		plt.clf()

		hhc_pitch_array = blade_results['hhc_pitch_array'][rotor_idx, blade_idx,:]
		
		#print(f'blade_twist_array: {blade_twist_array}')
		plt.figure(num=1)
		plt.plot(blade_twist_azimuth, hhc_pitch_array, linewidth=0.5)
		plt.title(f'{TITLE_DICT[plot_name.upper()]}: HHC over time.')
		plt.xlim([0, 360])
		plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name}_{ROTOR[rotor_idx]}_hhc.pdf', dpi=500, bbox_inches="tight", pad_inches=0.0)
		plt.cla()
		plt.clf()

		blade_flapping_array = blade_results['blade_flapping_array'][rotor_idx, blade_idx,:]
		
		#print(fblade_twist_array: {blade_twist_array}')
		plt.figure(num=1)
		plt.plot(blade_twist_azimuth, blade_flapping_array*(180.0/math.pi), linewidth=0.5)
		plt.title(f'{TITLE_DICT[plot_name.upper()]}: Flapping over time.')
		plt.xlim([0, 360])
		plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name}_{ROTOR[rotor_idx]}_flap.pdf', dpi=500, bbox_inches="tight", pad_inches=0.0)
		plt.cla()
		plt.clf()

		blade_flapping_der_array = blade_results['blade_flapping_der_array'][rotor_idx, blade_idx,:]
		
		#print(f'blade_twist_array: {blade_twist_array}')
		plt.figure(num=1)
		plt.plot(blade_twist_azimuth, blade_flapping_der_array, linewidth=0.5)
		plt.title(f'{TITLE_DICT[plot_name.upper()]}: Flapping derivative over time.')
		plt.xlim([0, 360])
		plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name}_{ROTOR[rotor_idx]}_flap_dot.pdf', dpi=500, bbox_inches="tight", pad_inches=0.0)
		plt.cla()
		plt.clf()

def plot_tr_wake(plot_name: str):
	blade_results = loadmat(f'{os.path.dirname(os.path.realpath(__file__))}/{plot_name.upper()}/results.mat')
	blade_results_hybrid = loadmat(f'{os.path.dirname(os.path.realpath(__file__))}/../helinovi_hybrid/helinovi/{plot_name.upper()}/results.mat')

	wake_filament_x = blade_results['wake_1_trajectory'][0,0,:]
	wake_filament_z = blade_results['wake_1_trajectory'][0,2,:]

	hybrid_wake_filament_x = blade_results_hybrid['wake_1_trajectory'][0,0,:]
	hybrid_wake_filament_z = blade_results_hybrid['wake_1_trajectory'][0,2,:]

	plt.figure(num = 1)

	plt.plot(wake_filament_x, wake_filament_z, 'b', hybrid_wake_filament_x, hybrid_wake_filament_z, 'k', linewidth=0.5)
	plt.ylim(-0.6, 1.25)
	plt.xlim(2.5, 5)
	plt.axis('scaled')
	plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name}_TR_wake_compare.pdf', dpi=500, bbox_inches="tight", pad_inches=0.0)
	plt.cla()
	plt.clf()
	

if __name__ == "__main__":

	matplotlib.rcParams['font.family'] = 'serif'
	matplotlib.rcParams['font.serif'] = prop.get_name()
	matplotlib.rcParams['mathtext.fontset'] = 'stix'

	#plot_blade_twist('BL')
	#plot_blade_twist('MN')
	
	#plot_mr_wake("ID2")
	#plot_mr_wake("ID1")

	# # # print("Plotting ID2_ASD")
	# # # #plot_blade_twist('ID2_ASD')
	# # # plot_blade_normal_pressures('ID2_ASD')
	# # # plot_acoustic_contours_fs("ID2_ASD")

	#print("Plotting ID1_ASU")
	#plot_hybrid_wake('ID1_ASU', 590)
	# # # # #plot_blade_twist('ID1_ASU')
	#plot_blade_normal_pressures('ID1_ASU', True)
	# # plot_acoustic_contours_fs("ID1_ASU")
	# # # # plot_acoustic_contours_mr("ID1_ASU")
	# plot_acoustic_contours_tr("ID1_ASU")
	#plot_acoustic_contours_tr("ID1_ASU", True)

	#print("Plotting ID1_ASD")
	#plot_hybrid_wake('ID1_ASD', 1050)
	# # # # #plot_blade_twist('ID1_ASD')
	#plot_blade_normal_pressures('ID1_ASD', True)
	# # # # #plot_wake_trajectory('BL')
	#plot_acoustic_contours_fs("ID1_ASD")
	# # # # plot_acoustic_contours_mr("ID1_ASD")
	# plot_acoustic_contours_tr("ID1_ASD")
	#plot_acoustic_contours_tr("ID1_ASD", True)

	# print("Plotting ID5")
	# # # # # # plot_blade_twist('ID5')
	# plot_blade_normal_pressures('ID5')
	# # # # # # #plot_wake_trajectory('MV')
	# plot_acoustic_contours_fs("ID5")
	# # # # # #plot_acoustic_contours_mr("ID5")
	# # # # # #plot_acoustic_contours_tr("ID5")

	#print("Plotting ID1")
	#plot_hybrid_wake('ID1', 173)
	# # # # plot_blade_twist('ID1')
	# plot_tr_wake('ID1')
	#plot_blade_normal_pressures('ID1', True, True)
	# # # # #plot_wake_trajectory('BL')
	#plot_acoustic_contours_fs("ID1", True)
	# plot_acoustic_contours_fs("ID1")
	
	# # # plot_acoustic_contours_mr("ID1")
	# plot_acoustic_contours_tr("ID1")

	print("Plotting ID2")
	#plot_hybrid_wake('ID2', 870)
	# # # plot_blade_twist('ID2')
	plot_blade_normal_pressures('ID2')
	plot_tr_wake('ID2')
	#plot_blade_normal_pressures('ID2', True, True)
	# #plot_wake_trajectory('BL')
	plot_acoustic_contours_fs("ID2")
	#plot_acoustic_contours_fs("ID2", True)
	# # plot_acoustic_contours_mr("ID2")
	# plot_acoustic_contours_tr("ID2")
	
	# plot_blade_twist('ID2_ASU')
	# plot_blade_normal_pressures('ID2_ASU')
	# #plot_wake_trajectory('MN')
	# plot_acoustic_contours_fs("ID2_ASU")
	# #plot_acoustic_contours_mr("ID2_ASU")
	# #plot_acoustic_contours_tr("ID2_ASU")


