#!/usr/bin/env python3
# Copyright (C) 2022, Robert F. Rau II, The Pennsylvania State University - All Rights Reserved
# You may use, distribute and modify this code under the
# terms of the GPLv3
#
# You should have received a copy of the GPLv3 license with
# this file. If not, please write to: rfr5313@psu.edu

import sys
import os

import subprocess
import argparse

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


def write_xfoil_inputs(Re: float, n: float, NACA: str | None, af_file: str | None, aoa_min: float, aoa_max: float, aoa_step: float):
	xfoil_input_filename = ""
	polar_output_filename = ""
	if NACA is not None:
		xfoil_input_filename = f"NACA{NACA}_{int(round(Re))}.xfi"
		polar_output_filename = f"NACA{NACA}_{int(round(Re))}_polar.dat"
	else:
		af_name = af_file.split("/")[-1].split(".")[0]
		xfoil_input_filename = f"{af_name}_{int(round(Re))}.xfi"
		polar_output_filename = f"{af_name}_{int(round(Re))}_polar.dat"

	with open(xfoil_input_filename, 'w') as xfoil_file:

		if NACA is not None:
			xfoil_file.write(f"NACA {NACA}\n")
		else:
			xfoil_file.write(f"load {af_file}\n")
		
		xfoil_file.write("pane\n")
		xfoil_file.write("oper\n")
		xfoil_file.write(f"visc {Re}\n")
		xfoil_file.write("iter 10000\n")
		xfoil_file.write("pacc \n\n\n")
		xfoil_file.write(f"aseq 0 {aoa_max} {aoa_step}\n")
		xfoil_file.write("init\n")
		xfoil_file.write(f"aseq 0 {aoa_min} {aoa_step}\n")
		xfoil_file.write("pacc\n")
		xfoil_file.write("psor 1\n")
		xfoil_file.write("pwrt 1\n")
		xfoil_file.write(f"{polar_output_filename}\n")
		xfoil_file.write("Y\n\n")
		xfoil_file.write("quit\n")

	return (xfoil_input_filename, polar_output_filename)

def run_xfoil(xfoil_path: str | None, xf_input_file: str):

	run_args = [f"{os.path.dirname(os.path.realpath(__file__))}/run_xfoil.bash"]

	if xfoil_path is not None:
		run_args.append(f"{xfoil_path}/")

	xf_output_file = xf_input_file.split(".")[0]+".out"

	run_args.append(xf_input_file)
	run_args.append(xf_output_file)

	subprocess.run(run_args)

def check_aerodas(polar_filename: str, tc_ratio: float, NACA: str | None, plot_suffix: str = "", plot_legend: bool = True):
	sys.path.insert(0, f'{os.path.dirname(os.path.realpath(__file__))}/../')
	import libopencopter as oc
	import matplotlib.pyplot as plt
	import numpy as np
	import math

	aerodas_model = oc.create_aerodas_from_xfoil_polar(polar_filename, tc_ratio)

	aoa_deg = np.linspace(-20, 110, 1000)
	aoa = (math.pi/180.0)*aoa_deg
	
	c_l = np.zeros(aoa.shape)
	c_d = np.zeros(aoa.shape)

	for idx, aoa in enumerate(aoa):
		c_l[idx] = aerodas_model.get_Cl(aoa, 0)
		c_d[idx] = aerodas_model.get_Cd(aoa, 0)

	xfoil_label = f'NACA {NACA} XFOIL data' if NACA is not None else "XFOIL data"

	print(f'plotting curves')
	plt.figure()
	plt.plot(aoa_deg, c_l, 'r', label = "AERODAS MODEL", linewidth=1)
	plt.plot(aerodas_model.alpha, aerodas_model.CL, 'b', label = xfoil_label, linewidth=1)
	plt.xlabel("$\\alpha$ [$^\circ$]", fontsize=font_size0)
	plt.ylabel("$C_l$", fontsize=font_size0)
	if plot_legend:
		plt.legend()
	plt.grid()
	plt.ylim(-1.5, 1.5)
	plt.xticks(fontsize=font_size3)
	plt.yticks(fontsize=font_size3)
	if plot_suffix != "":
		plt.savefig(f"Cl_{plot_suffix}.png", dpi=500, bbox_inches="tight", pad_inches=0.1)
	else:
		plt.savefig("Cl.png", dpi=500, bbox_inches="tight", pad_inches=0.1)

	plt.figure()
	plt.plot(aoa_deg, c_d, 'r', label = "AERODAS MODEL", linewidth=1)
	plt.plot(aerodas_model.alpha, aerodas_model.CD, 'b', label=xfoil_label, linewidth=1)
	plt.xlabel("$\\alpha$ [$^\circ$]", fontsize=font_size0)
	plt.ylabel("$C_d$", fontsize=font_size0)
	#plt.legend()
	plt.grid()
	plt.ylim(-0.1, 2.1)
	plt.xticks(fontsize=font_size3)
	plt.yticks(fontsize=font_size3)
	if plot_suffix != "":
		plt.savefig(f"Cd_{plot_suffix}.png", dpi=500, bbox_inches="tight", pad_inches=0.1)
	else:
		plt.savefig("Cd.png", dpi=500, bbox_inches="tight", pad_inches=0.1)
	#plt.show()

def generate_polar(Re, n, NACA, filename, aoa_min, aoa_max, xf):
	xfoil_input_file, polar_output_filename = write_xfoil_inputs(Re, n, NACA, filename, aoa_min, aoa_max)

	run_xfoil(xf, xfoil_input_file)
	return polar_output_filename


def main():

	matplotlib.rcParams['font.family'] = 'serif'
	matplotlib.rcParams['font.serif'] = prop.get_name()
	matplotlib.rcParams['mathtext.fontset'] = 'stix'

	parser = argparse.ArgumentParser("xf_gen", description="Xfoil data generator for AERODAS model input")

	parser.add_argument(
		"-f",
		type=str,
		help="Airfoil coordinate file, Selig format expected",
		required=False
	)

	parser.add_argument(
		"-N",
		type=str,
		help="NACA Airfoil number",
		required=False
	)

	parser.add_argument(
		"-re",
		type=float,
		help="Reynolds number to simulate airfoil at",
		required=True
	)

	parser.add_argument(
		"-n",
		type=float,
		help="Turbulent transition exponent, Ncrit. See Xfoil documentation for reasonable ranges. Default = 9",
		default=9.0,
		required=False
	)

	parser.add_argument(
		"-xf",
		type=str,
		help="Path to the Xfoil executable, if left blank, the executable is assumed to be in the PATH environment variable",
		required=False
	)

	parser.add_argument(
		"-aoa_min",
		type=float,
		help="Minimum angle of attack to perform sweep over. Default = -10 degrees",
		default=-10,
		required=False
	)

	parser.add_argument(
		"-aoa_max",
		type=float,
		help="Maximum angle of attack to perform sweep over. Default = 15 degrees",
		default=15,
		required=False
	)
	
	parser.add_argument(
		"-c",
		action="store_true",
		help="Check AeroDAS curve fitting from data",
		default=False,
		required=False
	)

	parser.add_argument(
		"-tc_ratio",
		type=float,
		help="Thickness to chord ratio of the current airfoil",
		required=False
	)

	parser.add_argument(
		"-d",
		action="store_true",
		help="dry run, don't actually generate polar. Use in conjunction with -c to check previously generated polar",
		default=False,
		required=False
	)
	
	parser.add_argument(
		"-aoa_step",
		type=float,
		help="Maximum angle of attack to perform sweep over. Default = 15 degrees",
		default=0.5,
		required=False
	)
	
	parser.add_argument(
		"-ps",
		type=str,
		help="Suffix for plot output",
		default="",
		required=False
	)

	parser.add_argument(
		"-l",
		action="store_true",
		help="plot legend",
		default=False,
		required=False
	)

	args = parser.parse_args()

	polar_output_filename = ""
	if args.N is not None:
		polar_output_filename = f"NACA{args.N}_{int(round(args.re))}_polar.dat"
	else:
		af_name = args.f.split("/")[-1].split(".")[0]
		polar_output_filename = f"{af_name}_{int(round(args.re))}_polar.dat"

	if not args.d:
		polar_output_filename = generate_polar(args.re, args.n, args.N, args.f, args.aoa_min, args.aoa_max, args.xf, args.aoa_step)

	if args.c:
		if args.tc_ratio is None:
			print("ERROR: The thickness to chord ratio of the airfoil must be specified to check the AERODAS model")
			quit()
		
		print("Checking aerodas model. Stand by for plots")
		check_aerodas(polar_output_filename, args.tc_ratio, args.N, args.ps, args.l)


if __name__ == "__main__":
	main()
