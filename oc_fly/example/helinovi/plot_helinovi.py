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
	"ID1": [(82, 100, 19), (82, 100, 19), (82, 100, 19)],
	"ID1_ASU": [(82, 100, 19), (82, 100, 19), (82, 100, 19)],
	"ID1_ASD": [(82, 100, 19), (82, 100, 19), (82, 100, 19)],
	"ID2": [(81, 99, 19), (81, 99, 19), (84, 101, 18)],
	"ID2_ASU": [(81, 99, 19), (81, 99, 19), (81, 99, 19)],
	"ID2_ASD": [(81, 99, 19), (81, 99, 19), (81, 99, 19)],
	"ID5": [(85, 119, 18), (85, 119, 18), (85, 114, 18)]
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
	"ID1": [(82, 100, 19), (82, 100, 19), (82, 100, 19)],
	"ID1_ASU": [(82, 100, 19), (82, 100, 19), (82, 100, 19)],
	"ID1_ASD": [(82, 100, 19), (82, 100, 19), (82, 100, 19)],
	"ID2": [(81, 99, 19), (81, 99, 19), (81, 99, 19)],
	"ID2_ASU": [(81, 99, 19), (81, 99, 19), (81, 99, 19)],
	"ID2_ASD": [(81, 99, 19), (81, 99, 19), (81, 99, 19)],
	"ID5": [(85, 119, 18), (85, 119, 18), (85, 119, 18)]
}

#cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["b", "dodgerblue", "cyan", "chartreuse", "yellow", "gold", "orange", "darkorange", "red"])
#cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["b", "dodgerblue", "cyan", "springgreen", "chartreuse", "greenyellow","yellow", "orange", "red"])
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["b", "dodgerblue", "cyan", "springgreen", "chartreuse", "yellow", "orange", "red"])
cmap_lines = matplotlib.colors.LinearSegmentedColormap.from_list("", ["black", "black"])

def plot_acoustic_contours_fs(plot_name: str):
	#x_grid, y_grid, measured = read_hart_contour_tecplot(f'{os.path.dirname(os.path.realpath(__file__))}/{plot_name.lower()}-contour-meas.tec')

	namelist = parse_namelist(f'{os.path.dirname(os.path.realpath(__file__))}/{plot_name.upper()}/acoustics/full_system/case.nam')
	wopwop_results = parse_wopwop_results_namelist(f'{os.path.dirname(os.path.realpath(__file__))}/{plot_name.upper()}/acoustics/full_system', namelist)

	i_max = len(wopwop_results.oaspl_dba_grid.obs_x)
	j_max = len(wopwop_results.oaspl_dba_grid.obs_x[0])

	#print(f'i_max: {i_max}')
	#print(f'j_max: {j_max}')

	phi = np.linspace(0, 2*math.pi, 1000)
	rotor_x = 2*np.cos(phi)
	rotor_z = 2*np.sin(phi)

	for r_idx, freq_range in enumerate(namelist.observers[0].ranges):
		oaspl_linear = [oaspl_db.functions[0].data[0] for oaspl_db in wopwop_results.frequency_ranges_db[r_idx]]

		oaspl_db = [[oaspl_linear[j*i_max + i] for i in range(i_max)] for j in range(j_max)]

		oaspl_linear = [oaspl_db.functions[0].data[0] for oaspl_db in wopwop_results.frequency_ranges_dba[r_idx]]

		oaspl_dba = [[oaspl_linear[j*i_max + i] for i in range(i_max)] for j in range(j_max)]

		#print(wopwop_results.oaspl_db_grid.obs_y[0])
		#print(wopwop_results.oaspl_db_grid.obs_x[1])

		y = wopwop_results.oaspl_dba_grid.obs_y[0]
		x = [_x[0] for _x in wopwop_results.oaspl_dba_grid.obs_x]
		x_delta = 4 - x[0]
		x.reverse()
		x = np.asarray(x)
		x = x + x_delta

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
		#ax0 = plt.subplot(121)
		#plt2 = plt.contour([-_y/R for _y in y], [-(_x - offset)/R for _x in x], oaspl_db, levels=clevels, linewidths=0.1)
		#plt1 = plt.contourf([-_y/R for _y in y], [-(_x - offset)/R for _x in x], oaspl_db, levels=clevels)
		#plt2 = plt.contour([-_y/R for _y in y], [-(_x - offset)/R for _x in x], oaspl_db, linewidths=0.1)
		#plt1 = plt.contourf([-_y/R for _y in y], [-(_x - offset)/R for _x in x], oaspl_db)
		plt2 = plt.contour(y, x, oaspl_db, levels=clevels, linewidths=0.5, cmap=cmap_lines, norm=cnorm)
		plt1 = plt.contourf(y, x, oaspl_db, levels=clevels, cmap=cmap, norm=cnorm)
		#plt.clabel(plt2, clevels, inline=True, colors='k', fontsize=5)
		clabels = plt.clabel(plt2, clevels, colors='k', fontsize=font_size4)
		plt.ylabel('x [m]')
		plt.xlabel('y [m]')
		plt.title('Prediction')
		plt.axis('scaled')
		#plt.axis('equal')

		for label in clabels:
			label.set_bbox(dict(facecolor='white', edgecolor='white', pad=0.01))

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
		plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name.upper()}_{freq_range.Title}_dB.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
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
		ax0 = plt.subplot(121)
		#plt2 = plt.contour([-_y/R for _y in y], [-(_x - offset)/R for _x in x], oaspl_db, levels=clevels, linewidths=0.1)
		#plt1 = plt.contourf([-_y/R for _y in y], [-(_x - offset)/R for _x in x], oaspl_db, levels=clevels)
		#plt2 = plt.contour([-_y/R for _y in y], [-(_x - offset)/R for _x in x], oaspl_db, linewidths=0.1)
		#plt1 = plt.contourf([-_y/R for _y in y], [-(_x - offset)/R for _x in x], oaspl_db)
		plt2 = plt.contour(y, x, oaspl_dba, levels=clevels, linewidths=0.5, cmap=cmap_lines, norm=cnorm)
		plt1 = plt.contourf(y, x, oaspl_dba, levels=clevels, cmap=cmap, norm=cnorm)
		plt.plot(rotor_x, rotor_z, "k", linewidth=1)
		#plt.clabel(plt2, clevels, inline=True, colors='k', fontsize=5)
		clabels = plt.clabel(plt2, clevels, colors='k', fontsize=font_size4)
		plt.ylabel('x [m]')
		plt.xlabel('y [m]')
		
		plt.title('Prediction')
		plt.axis('scaled')
		plt.ylim(-4, 4)
		#plt.axis('equal')

		for label in clabels:
			label.set_bbox(dict(facecolor='white', edgecolor='white', pad=0.01))

		ax1 = plt.subplot(122)
		arr=plt.imread('/mnt/c/Users/rfr2/Pictures/helinovi_level_flight.JPG')
		#arr=plt.imread('/mnt/c/Users/rfr2/Pictures/helinovi_level_flight_low_res-MikFaj7c_-transformed.png')
		
		plt.imshow(np.flipud(arr) ,interpolation='bilinear', origin='lower', extent=[y[0],y[-1],-4,4])
		plt.plot(rotor_x, rotor_z, "k", linewidth=1)
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
		plt.xlabel('y [m]')
		#plt.ylabel('-x/R')

		#ax = plt.subplot(133)
		cax = fig.add_axes([0.95, 0.2, 0.02, 0.582])
		fig.colorbar(plt1, cax, orientation='vertical', label='OASPL [dBA]')
		plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name.upper()}_{freq_range.Title}_dBA.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
		#plt.show()

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
		plt2 = plt.contour(y, x, oaspl_db, levels=15, linewidths=0.5, cmap=cmap, norm=cnorm)
		plt1 = plt.contourf(y, x, oaspl_db, levels=15, cmap=cmap, norm=cnorm)
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
		plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name.upper()}_{freq_range.Title}_MR_dB.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
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
		plt2 = plt.contour(y, x, oaspl_dba, levels=15, linewidths=0.1, cmap=cmap, norm=cnorm)
		plt1 = plt.contourf(y, x, oaspl_dba, levels=15, cmap=cmap, norm=cnorm)
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
		plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name.upper()}_{freq_range.Title}_MR_dBA.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
		#plt.show()

def plot_acoustic_contours_tr(plot_name: str):
	#x_grid, y_grid, measured = read_hart_contour_tecplot(f'{os.path.dirname(os.path.realpath(__file__))}/{plot_name.lower()}-contour-meas.tec')

	namelist = parse_namelist(f'{os.path.dirname(os.path.realpath(__file__))}/{plot_name.upper()}/acoustics/rotor_1/case.nam')
	wopwop_results = parse_wopwop_results_namelist(f'{os.path.dirname(os.path.realpath(__file__))}/{plot_name.upper()}/acoustics/rotor_1', namelist)

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
		plt2 = plt.contour(y, x, oaspl_db, levels=clevels, linewidths=0.1, cmap=cmap, norm=cnorm)
		plt1 = plt.contourf(y, x, oaspl_db, levels=clevels, cmap=cmap, norm=cnorm)
		#plt2 = plt.contour(y, x, oaspl_db, levels=15, linewidths=0.1)
		#plt1 = plt.contourf(y, x, oaspl_db, levels=15)	
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
		plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name.upper()}_{freq_range.Title}_TR_dB.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
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
		plt2 = plt.contour(y, x, oaspl_dba, levels=clevels, linewidths=0.1)
		plt1 = plt.contourf(y, x, oaspl_dba, levels=clevels)
		#plt2 = plt.contour(y, x, oaspl_dba, levels=15, linewidths=0.1)
		#plt1 = plt.contourf(y, x, oaspl_dba, levels=15)
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
		plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name.upper()}_{freq_range.Title}_TR_dBA.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
		#plt.show()

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

MEASURED_LOADING_TR = {
	"ID1": ["ID1_TR_loading_80_1r.tsv", "xx", "ID1_TR_loading_97_1r.tsv"],
	"ID1_ASU": ["xx", "xx", "ID1_ASU_TR_loading_97_1r.tsv"],
	"ID1_ASD": ["xx", "xx", "ID1_ASD_TR_loading_97_1r.tsv"],
	"ID2": ["ID2_TR_loading_80_1r.tsv", "xx", "ID2_TR_loading_97_1r.tsv"],
	"ID2_ASU": ["xx", "xx", "xx"],
	"ID2_ASD": ["xx", "xx", "xx"],
	"ID5": ["xx", "xx", "xx"]
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
	for s_idx in range(2):
		wake_x = wake_element_piv[0, s_idx, 0, :]
		wake_z = wake_element_piv[0, s_idx, 1, :]
		#print(f"Plotting slice {s_idx}")
		#print(wake_x[s_idx, :])
		#print(wake_z[s_idx, :])
		plt.figure(num = 1)
		plt.plot(rotor_x, rotor_z, 'k', -wake_x[:], wake_z[:], "r-", markersize=0.8)
		plt.axis('square')
		plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name}_{ROTOR[0]}_wake_trajectory_plane{s_idx}.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
		plt.cla()
		plt.clf()

def plot_blade_normal_pressures(plot_name: str):

	with open(f'{os.path.dirname(os.path.realpath(__file__))}/../helinovi_params.json5') as param_file:
		params = pyjson5.load(param_file)

	flight_condition = list(filter(lambda x: x["name"] == plot_name, params['flight_conditions']))[0]

	omegas = []
	for m in flight_condition["motion"]:
		if ('blade_element_func' in m) and (m['blade_element_func'] == "rotation"):
			omegas.append(m['omega'])
	
	blade_results = loadmat(f'{os.path.dirname(os.path.realpath(__file__))}/{plot_name.upper()}/results.mat')

	#measured_blade_loading = read_hart_blade_data_tecplot(f'{os.path.dirname(os.path.realpath(__file__))}/CnM2_{plot_name}_ca.tec')

	#measured_mean = measured_blade_loading.mean()
	#measured_blade_loading = measured_blade_loading - measured_mean
	
	for rotor_idx, rotor_blade_loading in enumerate(blade_results['span_element_af_loading']):
		for span_idx, blade_loading in enumerate(rotor_blade_loading):
			#blade_loading = blade_results['span_element_loading'][0][0]

			omega = omegas[rotor_idx]

			#blade_loading = blade_loading - computed_mean

			measured_blade_loading = np.asarray([])#np.zeros(1440)
			measured_azimuth = np.asarray([])#np.linspace(0, 360, measured_blade_loading.size)

			if os.path.exists(f'{os.path.dirname(os.path.realpath(__file__))}/{MEASURED_LOADING[plot_name][rotor_idx][span_idx]}'):
				measured_loading = read_hart_blade_data_tsv(f'{os.path.dirname(os.path.realpath(__file__))}/{MEASURED_LOADING[plot_name][rotor_idx][span_idx]}')
				measured_azimuth = measured_loading[0,:]
				measured_blade_loading = measured_loading[1,:]
				#measured_blade_loading = measured_blade_loading - measured_blade_loading.mean()

			#computational_azimuth = np.linspace(0, 360, blade_loading.size)
			computational_azimuth = blade_results['span_element_azimuth'][rotor_idx]

			print(f'computational_azimuth: {computational_azimuth}')
			#measured_blade_loading = np.zeros(measured_azimuth.size)

			delta_azi = 0
			# if (rotor_idx == 1) and (plot_name != "ID1_ASU"):
			# 	delta_azi = 25
			# elif(rotor_idx == 1) and (plot_name == "ID1_ASU"):
			# 	delta_azi = 180

			plt.figure(num=1)
			plt.plot(measured_azimuth, measured_blade_loading - measured_blade_loading.mean(), 'ro-', computational_azimuth + delta_azi, blade_loading - blade_loading.mean(), 'b-', linewidth=0.5, markersize=1.0, mfc='none')
			plt.title(f'{TITLE_DICT[plot_name.upper()]} {SPAN_LOCATION[span_idx]*100}% Span blade loading. Mean removed.')
			plt.legend(['Measured', 'OpenCOPTER'])
			plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name}_{ROTOR[rotor_idx]}_{SPAN_LOCATION[span_idx]:0.2f}_normal_loading.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
			plt.cla()
			plt.clf()

			if rotor_idx == 1:
				start_azi = 360 - delta_azi
				measured_azimuth_tr = []
				measured_blade_loading_tr = []
				if os.path.exists(f'{os.path.dirname(os.path.realpath(__file__))}/{MEASURED_LOADING_TR[plot_name][span_idx]}'):
					measured_loading_tr = read_hart_blade_data_tsv(f'{os.path.dirname(os.path.realpath(__file__))}/{MEASURED_LOADING_TR[plot_name][span_idx]}')
					measured_azimuth_tr = measured_loading_tr[0,:]
					measured_blade_loading_tr = measured_loading_tr[1,:]

				plt.figure(num=1)
				plt.plot(measured_azimuth_tr, measured_blade_loading_tr - np.asarray(measured_blade_loading_tr).mean(), computational_azimuth[0:360], blade_loading[start_azi:start_azi+360] - np.asarray(blade_loading[start_azi:start_azi+360]).mean(), linewidth=0.5)
				plt.title(f'{TITLE_DICT[plot_name.upper()]} {SPAN_LOCATION[span_idx]*100}% Span blade loading. Mean removed.')
				plt.legend(['Measured', 'OpenCOPTER'])
				plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name}_{ROTOR[rotor_idx]}_{SPAN_LOCATION[span_idx]:0.2f}_normal_loading_1rev.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
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

			plt.figure(num=1)
			plt.plot(measured_azimuth, measured_blade_loading, computational_azimuth, blade_loading, linewidth=0.5)
			plt.title(f'{TITLE_DICT[plot_name.upper()]}: 87% Span blade loading. 10/rev highpassed.')
			plt.legend(['Measured', 'OpenCOPTER'])
			#plt.xlim([20, 90])
			plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name}_{ROTOR[rotor_idx]}_{SPAN_LOCATION[span_idx]:0.2f}_normal_loading_hp.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
			plt.cla()
			plt.clf()

			aoa_eff = blade_results['span_element_aoa_eff'][rotor_idx][span_idx]
			
			plt.figure(num=1)
			plt.plot(computational_azimuth, aoa_eff, linewidth=0.5)
			plt.title(f'{TITLE_DICT[plot_name.upper()]} 87% Effective AOA.')
			#plt.legend(['Measured', 'OpenCOPTER'])
			plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name}_{ROTOR[rotor_idx]}_{SPAN_LOCATION[span_idx]:0.2f}_aoa_eff.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
			plt.cla()
			plt.clf()

			u_p = blade_results['span_element_up'][rotor_idx][span_idx]
			
			plt.figure(num=1)
			plt.plot(computational_azimuth, u_p, linewidth=0.5)
			plt.title(f'{TITLE_DICT[plot_name.upper()]} 87% U_p.')
			#plt.legend(['Measured', 'OpenCOPTER'])
			plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name}_{ROTOR[rotor_idx]}_{SPAN_LOCATION[span_idx]:0.2f}_u_p.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
			plt.cla()
			plt.clf()

			aoa = blade_results['span_element_aoa'][rotor_idx][span_idx]
			
			plt.figure(num=1)
			plt.plot(computational_azimuth, aoa, linewidth=0.5)
			plt.title(f'{TITLE_DICT[plot_name.upper()]} 87% AOA.')
			#plt.legend(['Measured', 'OpenCOPTER'])
			plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name}_{ROTOR[rotor_idx]}_{SPAN_LOCATION[span_idx]:0.2f}_aoa.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
			plt.cla()
			plt.clf()

			theta = blade_results['span_element_theta'][rotor_idx][span_idx]
			
			plt.figure(num=1)
			plt.plot(computational_azimuth, theta, linewidth=0.5)
			plt.title(f'{TITLE_DICT[plot_name.upper()]} 87% theta.')
			#plt.legend(['Measured', 'OpenCOPTER'])
			plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name}_{ROTOR[rotor_idx]}_{SPAN_LOCATION[span_idx]:0.2f}_theta.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
			plt.cla()
			plt.clf()

			inflow_angle = blade_results['span_element_inflow_angle'][rotor_idx][span_idx]
			
			plt.figure(num=1)
			plt.plot(computational_azimuth, inflow_angle, linewidth=0.5)
			plt.title(f'{TITLE_DICT[plot_name.upper()]} 87% inflow angle.')
			#plt.legend(['Measured', 'OpenCOPTER'])
			plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name}_{ROTOR[rotor_idx]}_{SPAN_LOCATION[span_idx]:0.2f}_inflow_angle.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
			plt.cla()
			plt.clf()

			gamma = blade_results['span_element_gamma'][rotor_idx][span_idx]
			
			plt.figure(num=1)
			plt.plot(computational_azimuth, gamma, linewidth=0.5)
			plt.title(f'{TITLE_DICT[plot_name.upper()]} 87% gamma.')
			#plt.legend(['Measured', 'OpenCOPTER'])
			plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name}_{ROTOR[rotor_idx]}_{SPAN_LOCATION[span_idx]:0.2f}_gamma.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
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
		plt.plot(measured_x, measured_z, '*', x_tppc, z_tppc, '-', linewidth=0.5)
		plt.title(f'{TITLE_DICT[plot_name.upper()]}: Y slice: {target_y_slice:.2f}')
		plt.legend(['Measured', 'OpenCOPTER'])
		#plt.axis('scaled')
		plt.xlim([-1.1, 1.1])
		plt.ylim([-0.1, 0.15])
		plt.gca().invert_xaxis()
		plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/Hart_{plot_name}_wake_y{target_y_slice:.2f}.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
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
		plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/Hart_{plot_name}_core_size_y{target_y_slice:.2f}.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
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
		plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name}_{ROTOR[rotor_idx]}_blade_twist.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
		plt.cla()
		plt.clf()
		
		elastic_twist_array = blade_results['elastic_twist_array'][rotor_idx, blade_idx,:]
		
		#print(f'blade_twist_array: {blade_twist_array}')
		plt.figure(num=1)
		plt.plot(blade_twist_azimuth, elastic_twist_array, linewidth=0.5)
		plt.title(f'{TITLE_DICT[plot_name.upper()]}: Elastic twist over time.')
		plt.xlim([0, 360])
		plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name}_{ROTOR[rotor_idx]}_elastic_twist.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
		plt.cla()
		plt.clf()

		collective_pitch_array = blade_results['collective_pitch_array'][rotor_idx, :]
		
		#print(f'blade_twist_array: {blade_twist_array}')
		plt.figure(num=1)
		plt.plot(blade_twist_azimuth, collective_pitch_array, linewidth=0.5)
		plt.title(f'{TITLE_DICT[plot_name.upper()]}: Collective pitch over time.')
		plt.xlim([0, 360])
		plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name}_{ROTOR[rotor_idx]}_collective.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
		plt.cla()
		plt.clf()

		sin_pitch_array = blade_results['sin_pitch_array'][rotor_idx, blade_idx,:]
		sin_pitch_array_hp = sig.sosfilt(computed_sos, sin_pitch_array)

		#print(f'blade_twist_array: {blade_twist_array}')
		plt.figure(num=1)
		plt.plot(blade_twist_azimuth, sin_pitch_array, linewidth=0.5)
		plt.title(f'{TITLE_DICT[plot_name.upper()]}: 1S over time.')
		plt.xlim([0, 360])
		plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name}_{ROTOR[rotor_idx]}_1s.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
		plt.cla()
		plt.clf()

		#print(f'blade_twist_array: {blade_twist_array}')
		plt.figure(num=1)
		plt.plot(blade_twist_azimuth, sin_pitch_array_hp, linewidth=0.5)
		plt.title(f'{TITLE_DICT[plot_name.upper()]}: 1S over time HP.')
		plt.xlim([0, 360])
		plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name}_{ROTOR[rotor_idx]}_1s_hp.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
		plt.cla()
		plt.clf()

		cos_pitch_array = blade_results['cos_pitch_array'][rotor_idx, blade_idx,:]
		cos_pitch_array_hp = sig.sosfilt(computed_sos, cos_pitch_array)

		#print(f'blade_twist_array: {blade_twist_array}')
		plt.figure(num=1)
		plt.plot(blade_twist_azimuth, cos_pitch_array, linewidth=0.5)
		plt.title(f'{TITLE_DICT[plot_name.upper()]}: 1C over time.')
		plt.xlim([0, 360])
		plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name}_{ROTOR[rotor_idx]}_1c.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
		plt.cla()
		plt.clf()

		plt.figure(num=1)
		plt.plot(blade_twist_azimuth, cos_pitch_array_hp, linewidth=0.5)
		plt.title(f'{TITLE_DICT[plot_name.upper()]}: 1C over time HP.')
		plt.xlim([0, 360])
		plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name}_{ROTOR[rotor_idx]}_1c_hp.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
		plt.cla()
		plt.clf()

		hhc_pitch_array = blade_results['hhc_pitch_array'][rotor_idx, blade_idx,:]
		
		#print(f'blade_twist_array: {blade_twist_array}')
		plt.figure(num=1)
		plt.plot(blade_twist_azimuth, hhc_pitch_array, linewidth=0.5)
		plt.title(f'{TITLE_DICT[plot_name.upper()]}: HHC over time.')
		plt.xlim([0, 360])
		plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name}_{ROTOR[rotor_idx]}_hhc.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
		plt.cla()
		plt.clf()

		blade_flapping_array = blade_results['blade_flapping_array'][rotor_idx, blade_idx,:]
		
		#print(fblade_twist_array: {blade_twist_array}')
		plt.figure(num=1)
		plt.plot(blade_twist_azimuth, blade_flapping_array*(180.0/math.pi), linewidth=0.5)
		plt.title(f'{TITLE_DICT[plot_name.upper()]}: Flapping over time.')
		plt.xlim([0, 360])
		plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name}_{ROTOR[rotor_idx]}_flap.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
		plt.cla()
		plt.clf()

		blade_flapping_der_array = blade_results['blade_flapping_der_array'][rotor_idx, blade_idx,:]
		
		#print(f'blade_twist_array: {blade_twist_array}')
		plt.figure(num=1)
		plt.plot(blade_twist_azimuth, blade_flapping_der_array, linewidth=0.5)
		plt.title(f'{TITLE_DICT[plot_name.upper()]}: Flapping derivative over time.')
		plt.xlim([0, 360])
		plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/HeliNOVI_{plot_name}_{ROTOR[rotor_idx]}_flap_dot.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
		plt.cla()
		plt.clf()


if __name__ == "__main__":

	matplotlib.rcParams['font.family'] = 'serif'
	matplotlib.rcParams['font.serif'] = prop.get_name()
	matplotlib.rcParams['mathtext.fontset'] = 'stix'

	#plot_blade_twist('BL')
	#plot_blade_twist('MN')
	
	plot_mr_wake("ID2")

	#plot_blade_twist('ID2_ASD')
	plot_blade_normal_pressures('ID2_ASD')
	plot_acoustic_contours_fs("ID2_ASD")

	#plot_blade_twist('ID1_ASU')
	plot_blade_normal_pressures('ID1_ASU')
	plot_acoustic_contours_fs("ID1_ASU")
	plot_acoustic_contours_mr("ID1_ASU")
	plot_acoustic_contours_tr("ID1_ASU")

	#plot_blade_twist('ID1_ASD')
	plot_blade_normal_pressures('ID1_ASD')
	#plot_wake_trajectory('BL')
	plot_acoustic_contours_fs("ID1_ASD")
	plot_acoustic_contours_mr("ID1_ASD")
	plot_acoustic_contours_tr("ID1_ASD")

	# plot_blade_twist('ID5')
	# plot_blade_normal_pressures('ID5')
	# #plot_wake_trajectory('MV')
	# plot_acoustic_contours_fs("ID5")
	# #plot_acoustic_contours_mr("ID5")
	# #plot_acoustic_contours_tr("ID5")

	plot_blade_twist('ID1')
	plot_blade_normal_pressures('ID1')
	#plot_wake_trajectory('BL')
	plot_acoustic_contours_fs("ID1")
	plot_acoustic_contours_mr("ID1")
	plot_acoustic_contours_tr("ID1")

	# plot_blade_twist('ID2')
	plot_blade_normal_pressures('ID2')
	#plot_wake_trajectory('BL')
	plot_acoustic_contours_fs("ID2")
	# plot_acoustic_contours_mr("ID2")
	# plot_acoustic_contours_tr("ID2")
	
	# plot_blade_twist('ID2_ASU')
	# plot_blade_normal_pressures('ID2_ASU')
	# #plot_wake_trajectory('MN')
	# plot_acoustic_contours_fs("ID2_ASU")
	# #plot_acoustic_contours_mr("ID2_ASU")
	# #plot_acoustic_contours_tr("ID2_ASU")


