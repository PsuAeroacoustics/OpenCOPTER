#!/usr/bin/env python3

import sys
import os

#sys.path.insert(0, f'{os.path.dirname(os.path.realpath(__file__))}/../../../')
sys.path.insert(0, f'{os.path.dirname(os.path.realpath(__file__))}/../../dependencies/wopwopd')

#from libopencopter import *
from libwopwopd import *
import numpy as np
import math
import time
import argparse

import matplotlib
import matplotlib.pyplot as plt

R = 2

def read_hart_tecplot(filename, i_max = 17, j_max = 13):
	with open(filename, 'r') as tecfile:
		_ = tecfile.readline()
		_ = tecfile.readline()
		_ = tecfile.readline()

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

def plot_mv():
    x_grid, y_grid, measured = read_hart_tecplot(f'{os.path.dirname(os.path.realpath(__file__))}/mv-contour-meas.tec', 16, 13)

    wopwop_results = parse_wopwop_results(f'{os.path.dirname(os.path.realpath(__file__))}/MV/acoustics/full_system', 'case.nam')

    i_max = len(wopwop_results.oaspl_db_grid.obs_x)
    j_max = len(wopwop_results.oaspl_db_grid.obs_x[0])

    print(f'i_max: {i_max}')
    print(f'j_max: {j_max}')

    oaspl_linear = [oaspl_db.functions[1].data[0] for oaspl_db in wopwop_results.oaspl_db]

    oaspl_db = [[oaspl_linear[i*j_max + j] for i in range(i_max)] for j in range(j_max)]

    #print(wopwop_results.oaspl_db_grid.obs_y[0])
    #print(wopwop_results.oaspl_db_grid.obs_x[1])

    x = wopwop_results.oaspl_db_grid.obs_x[0]
    y = [_y[0] for _y in wopwop_results.oaspl_db_grid.obs_y]
    y.reverse()

    #print(x)
    offset = x[0] - -4.0

    clevels = np.linspace(85, 117, 17)

    #light_rainbow = cmap_map(lambda x: x/2 + 0.5, matplotlib.cm.rainbow)

    print([(_x - offset)/R for _x in x])
    print([_y/R for _y in y])
    print(x_grid[0,:])
    print(y_grid[:,0])
    #print(x)

    fig = plt.figure()
    ax0 = plt.subplot(121)
    plt2 = plt.contour([_y/R for _y in y], [(_x - offset)/R for _x in x], oaspl_db, levels=clevels, linewidths=0.1)
    plt1 = plt.contourf([_y/R for _y in y], [(_x - offset)/R for _x in x], oaspl_db, levels=clevels)
    #plt.clabel(plt2, clevels, inline=True, colors='k', fontsize=5)
    plt.clabel(plt2, clevels, colors='k', fontsize=5)
    plt.ylabel('x/R')
    plt.xlabel('y/R')
    plt.title('Prediction')
    #plt.axis('scaled')
    #plt.axis('equal')

    ax1 = plt.subplot(122)
    ax1.set_yticklabels([])
    #ax.set_xticklabels([])
    plt2 = plt.contour(x_grid[0,:], y_grid[:,0], measured, levels=clevels, linewidths=0.1)
    plt.contourf(x_grid[0,:], y_grid[:,0], measured, levels=clevels)
    plt.clabel(plt2, clevels, colors='k', fontsize=5)
    plt.title('Measured')
    #plt.clabel(CS, clevels, inline=True)
    #plt.contourf(y_grid[:,0], x_grid[0,:], measured, levels=clevels)
    #plt.axis('equal')
    #plt.colorbar()
    plt.xlabel('y/R')
    #plt.ylabel('-x/R')

    #ax = plt.subplot(133)
    cax = fig.add_axes([0.95, 0.11, 0.02, 0.77])
    fig.colorbar(plt1, cax, orientation='vertical', label='BVI SPL [dB]')
    plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/Hart_MV.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
    #plt.show()

def plot_mn():
    x_grid, y_grid, measured = read_hart_tecplot(f'{os.path.dirname(os.path.realpath(__file__))}/mn-contour-meas.tec')

    wopwop_results = parse_wopwop_results(f'{os.path.dirname(os.path.realpath(__file__))}/MN/acoustics/full_system', 'case.nam')

    i_max = len(wopwop_results.oaspl_db_grid.obs_x)
    j_max = len(wopwop_results.oaspl_db_grid.obs_x[0])

    print(f'i_max: {i_max}')
    print(f'j_max: {j_max}')

    oaspl_linear = [oaspl_db.functions[1].data[0] for oaspl_db in wopwop_results.oaspl_db]

    oaspl_db = [[oaspl_linear[i*j_max + j] for i in range(i_max)] for j in range(j_max)]

    #print(wopwop_results.oaspl_db_grid.obs_y[0])
    #print(wopwop_results.oaspl_db_grid.obs_x[1])

    x = wopwop_results.oaspl_db_grid.obs_x[0]
    y = [_y[0] for _y in wopwop_results.oaspl_db_grid.obs_y]
    y.reverse()

    #print(x)
    offset = x[0] - -4.0

    clevels = np.linspace(85, 117, 17)

    #light_rainbow = cmap_map(lambda x: x/2 + 0.5, matplotlib.cm.rainbow)

    print([(_x - offset)/R for _x in x])
    print([_y/R for _y in y])
    print(x_grid[0,:])
    print(y_grid[:,0])
    #print(x)

    fig = plt.figure()
    ax0 = plt.subplot(121)
    plt2 = plt.contour([_y/R for _y in y], [(_x - offset)/R for _x in x], oaspl_db, levels=clevels, linewidths=0.1)
    plt1 = plt.contourf([_y/R for _y in y], [(_x - offset)/R for _x in x], oaspl_db, levels=clevels)
    #plt.clabel(plt2, clevels, inline=True, colors='k', fontsize=5)
    plt.clabel(plt2, clevels, colors='k', fontsize=5)
    plt.ylabel('x/R')
    plt.xlabel('y/R')
    plt.title('Prediction')
    #plt.axis('scaled')
    #plt.axis('equal')

    ax1 = plt.subplot(122)
    ax1.set_yticklabels([])
    #ax.set_xticklabels([])
    plt2 = plt.contour(x_grid[0,:], y_grid[:,0], measured, levels=clevels, linewidths=0.1)
    plt.contourf(x_grid[0,:], y_grid[:,0], measured, levels=clevels)
    plt.clabel(plt2, clevels, colors='k', fontsize=5)
    plt.title('Measured')
    #plt.clabel(CS, clevels, inline=True)
    #plt.contourf(y_grid[:,0], x_grid[0,:], measured, levels=clevels)
    #plt.axis('equal')
    #plt.colorbar()
    plt.xlabel('y/R')
    #plt.ylabel('-x/R')

    #ax = plt.subplot(133)
    cax = fig.add_axes([0.95, 0.11, 0.02, 0.77])
    fig.colorbar(plt1, cax, orientation='vertical', label='BVI SPL [dB]')
    plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/Hart_MN.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
    #plt.show()

def plot_baseline():
    x_grid, y_grid, measured = read_hart_tecplot(f'{os.path.dirname(os.path.realpath(__file__))}/bl-contour-meas.tec')

    wopwop_results = parse_wopwop_results(f'{os.path.dirname(os.path.realpath(__file__))}/BL/acoustics/full_system', 'case.nam')

    i_max = len(wopwop_results.oaspl_db_grid.obs_x)
    j_max = len(wopwop_results.oaspl_db_grid.obs_x[0])

    print(f'i_max: {i_max}')
    print(f'j_max: {j_max}')

    oaspl_linear = [oaspl_db.functions[1].data[0] for oaspl_db in wopwop_results.oaspl_db]

    oaspl_db = [[oaspl_linear[i*j_max + j] for i in range(i_max)] for j in range(j_max)]

    #print(wopwop_results.oaspl_db_grid.obs_y[0])
    #print(wopwop_results.oaspl_db_grid.obs_x[1])

    x = wopwop_results.oaspl_db_grid.obs_x[0]
    y = [_y[0] for _y in wopwop_results.oaspl_db_grid.obs_y]
    y.reverse()

    #print(x)
    offset = x[0] - -4.0

    clevels = np.linspace(85, 117, 17)

    #light_rainbow = cmap_map(lambda x: x/2 + 0.5, matplotlib.cm.rainbow)

    print([(_x - offset)/R for _x in x])
    print([_y/R for _y in y])
    print(x_grid[0,:])
    print(y_grid[:,0])
    #print(x)

    fig = plt.figure()
    ax0 = plt.subplot(121)
    plt2 = plt.contour([_y/R for _y in y], [(_x - offset)/R for _x in x], oaspl_db, levels=clevels, linewidths=0.1)
    plt1 = plt.contourf([_y/R for _y in y], [(_x - offset)/R for _x in x], oaspl_db, levels=clevels)
    #plt.clabel(plt2, clevels, inline=True, colors='k', fontsize=5)
    plt.clabel(plt2, clevels, colors='k', fontsize=5)
    plt.ylabel('x/R')
    plt.xlabel('y/R')
    plt.title('Prediction')
    #plt.axis('scaled')
    #plt.axis('equal')

    ax1 = plt.subplot(122)
    ax1.set_yticklabels([])
    #ax.set_xticklabels([])
    plt2 = plt.contour(x_grid[0,:], y_grid[:,0], measured, levels=clevels, linewidths=0.1)
    plt.contourf(x_grid[0,:], y_grid[:,0], measured, levels=clevels)
    plt.clabel(plt2, clevels, colors='k', fontsize=5)
    plt.title('Measured')
    #plt.clabel(CS, clevels, inline=True)
    #plt.contourf(y_grid[:,0], x_grid[0,:], measured, levels=clevels)
    #plt.axis('equal')
    #plt.colorbar()
    plt.xlabel('y/R')
    #plt.ylabel('-x/R')

    #ax = plt.subplot(133)
    cax = fig.add_axes([0.95, 0.11, 0.02, 0.77])
    fig.colorbar(plt1, cax, orientation='vertical', label='BVI SPL [dB]')
    plt.savefig(f'{os.path.dirname(os.path.realpath(__file__))}/Hart_BL.png', dpi=500, bbox_inches="tight", pad_inches=0.0)
    #plt.show()

if __name__ == "__main__":
    plot_baseline()
    plot_mn()
    plot_mv()