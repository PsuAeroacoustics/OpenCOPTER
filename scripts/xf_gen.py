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

def write_xfoil_inputs(Re: float, n: float, NACA: str | None, af_file: str | None, aoa_min: float, aoa_max: float):
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
		xfoil_file.write(f"aseq 0 {aoa_max} 0.5\n")
		xfoil_file.write("init\n")
		xfoil_file.write(f"aseq 0 {aoa_min} 0.5\n")
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

def check_aerodas(polar_filename: str, tc_ratio: float):
	sys.path.insert(0, f'{os.path.dirname(os.path.realpath(__file__))}/../')
	import libopencopter as oc
	import matplotlib.pyplot as plt
	import numpy as np
	import math

	aerodas_model = oc.create_aerodas_from_xfoil_polar(polar_filename, tc_ratio)

	aoa_deg = np.linspace(-20, 110, 130)
	aoa = (math.pi/180.0)*aoa_deg
	
	c_l = np.zeros(aoa.shape)
	c_d = np.zeros(aoa.shape)

	for idx, aoa in enumerate(aoa):
		c_l[idx] = aerodas_model.get_Cl(aoa, 0)
		c_d[idx] = aerodas_model.get_Cd(aoa, 0)

	plt.figure()
	plt.plot(aoa_deg, c_l, label = "$AERODAS MODEL$")
	plt.plot(aerodas_model.alpha, aerodas_model.CL,label = "XFOIL data")
	plt.xlabel("Angle of attack (degrees)")
	plt.ylabel("Coefficient of lift")
	plt.legend()
	plt.grid()
	#plt.savefig("Cl.png", ["dpi": 500])

	plt.figure()
	plt.plot(aoa_deg, c_d, label = "$AERODAS MODEL$")
	plt.plot(aerodas_model.alpha, aerodas_model.CD, label="XFOIL data")
	plt.xlabel("Angle of attack (degrees)")
	plt.ylabel("Coefficient of drag")
	plt.legend()
	plt.grid()
	#plt.savefig("Cd.png", ["dpi": 500])
	plt.show()

def main():
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
	
	args = parser.parse_args()

	xfoil_input_file, polar_output_filename = write_xfoil_inputs(args.re, args.n, args.N, args.f, args.aoa_min, args.aoa_max)

	print("Running Xfoil")
	run_xfoil(args.xf, xfoil_input_file)

	if args.c:
		if args.tc_ratio is None:
			print("ERROR: The thickness to chord ratio of the airfoil must be specified to check the AERODAS model")
			quit()
		
		print("Checking aerodas model. Stand by for plots")
		check_aerodas(polar_output_filename, args.tc_ratio)


if __name__ == "__main__":
	main()
