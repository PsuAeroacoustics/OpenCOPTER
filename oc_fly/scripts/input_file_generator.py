#!/usr/bin/env python3
import json
import numpy as np
import math
import argparse
import copy

def main():

	parser = argparse.ArgumentParser("Input File Generator", description="Generate input files for the computational pipeline")

	parser.add_argument(
		"--base_geom",
		type=str,
		help="Base geometry file to work from",
		default=None,
		required=True
	)
	parser.add_argument(
		"--base_param",
		type=str,
		help="Base parameter file to work from",
		default=None,
		required=True
	)
	parser.add_argument(
		"--aoas",
		type=float,
		nargs="+",
		help="",
		default=[-5, 5, 2],
		required=False
	)
	parser.add_argument(
		"--aoa_direct",
		help="take h_sep as direct values, not a range spec",
		action='store_true',
		default=False,
		required=False
	)
	parser.add_argument(
		"--h_sep",
		type=float,
		nargs="+",
		help="",
		default=[-0.5, 12, 2],
		required=False
	)
	parser.add_argument(
		"--h_direct",
		help="take h_sep as direct values, not a range spec",
		action='store_true',
		default=False,
		required=False
	)
	parser.add_argument(
		"--v_sep",
		type=float,
		nargs="+",
		help="",
		default=[1, -2, 2],
		required=False
	)
	parser.add_argument(
		"--v_direct",
		help="take h_sep as direct values, not a range spec",
		action='store_true',
		default=False,
		required=False
	)
	parser.add_argument(
		"--l_sep",
		type=float,
		nargs="+",
		help="",
		default=[-3, 3, 2],
		required=False
	)
	parser.add_argument(
		"--l_direct",
		help="take h_sep as direct values, not a range spec",
		action='store_true',
		default=False,
		required=False
	)
	parser.add_argument(
		"--v_infs",
		type=float,
		nargs="+",
		help="",
		default=[30, 50, 2],
		required=False
	)
	parser.add_argument(
		"--v_inf_direct",
		help="take h_sep as direct values, not a range spec",
		action='store_true',
		default=False,
		required=False
	)

	parser.add_argument(
		"--c_t",
		type=float,
		nargs="+",
		help="",
		default=[0.005, 0.005, 1],
		required=False
	)
	parser.add_argument(
		"--c_t_direct",
		help="take c_t as direct values, not a range spec",
		action='store_true',
		default=False,
		required=False
	)


	args = parser.parse_args()

	if args.aoa_direct:
		aoas = np.asarray(args.aoas)
	else:
		aoas = np.linspace(args.aoas[0], args.aoas[1], int(args.aoas[2]))

	if args.h_direct:
		h_sep = np.asarray(args.h_sep)
	else:
		h_sep = np.geomspace(args.h_sep[0] + 1.0, args.h_sep[1] + 1.0, int(args.h_sep[2])) - 1.0

	if args.v_direct:
		v_sep = np.asarray(args.v_sep)
	else:
		v_sep = np.linspace(args.v_sep[0], args.v_sep[1], int(args.v_sep[2]))

	if args.l_direct:
		l_sep = np.asarray(args.l_sep)
	else:
		l_sep = np.linspace(args.l_sep[0], args.l_sep[1], int(args.l_sep[2]))

	if args.v_inf_direct:
		V_infs = np.asarray(args.v_infs)
	else:
		V_infs = np.linspace(args.v_infs[0], args.v_infs[1], int(args.v_infs[2]))

	if args.c_t_direct:
		c_ts = np.asarray(args.c_t)
	else:
		c_ts = np.linspace(args.c_t[0], args.c_t[1], int(args.c_t[2]))

	with open(args.base_geom) as geom_file:
		base_geom = json.load(geom_file)

	with open(args.base_param) as param_file:
		base_param = json.load(param_file)

	#num_rotors = 2
	#num_blades = 4
	#radius = 2
	#AR = 16.5
	#theta_tw_1 = -8.0
	#collective = 5.06*np.ones(num_rotors)
	#collective = collective.tolist()
	#aoa = np.linspace(-5,5,3).tolist()
	#aoa_rotors = [2,2]

	#omega = 109.12

	#V_inf = np.linspace(30.0,35.0,5).tolist()
	#h_sep = np.linspace(-2,2,5)
	#v_sep = np.linspace(-2,2,5)
	#l_sep = np.linspace(-2,2,5)

	print(base_geom)

	

	origin_idx = 1
	aircrafts = []
	a=0
	b=0
	c=0
	origin = np.zeros((h_sep.size*v_sep.size*l_sep.size, 3))
	#origin_aft_rotor = np.zeros((h_sep.size*v_sep.size*l_sep.size, 3))
	#origin_front_rotor = [0, 0, 0]
	for i in range(h_sep.size*v_sep.size*l_sep.size - 1):
		#print("a =",a," b =",b , " c =",c)
		origin_aft_rotor = [-(2 + h_sep[a]), l_sep[c], v_sep[b]]
		
		if(a==h_sep.size-1 and b!=v_sep.size-1):
			a = 0
			b = b+1
		elif(b==v_sep.size-1 and a==h_sep.size-1):
			a = 0
			b = 0
			c = c+1
		else:
			a =a+1
		origin[i] = origin_aft_rotor
	
		ac = copy.deepcopy(base_geom[0])
		ac["origin"][origin_idx] = origin_aft_rotor
		ac["Aircraft"] = f'{ac["Aircraft"]}_h{h_sep[a]:.2f}_v{v_sep[b]}_l{l_sep[c]}'

		aircrafts.append(ac)

	#origin = origin.tolist()


	# Data to be written in geometry file
	# Geom = {
	# 	"Aircraft": [
	# 		"Hart_II"
	# 	],
	# 	"AR": AR,
	# 	"collective": collective,
	# 	"number of blades": num_rotors,
	# 	"number of rotors": num_blades,
	# 	"origin": origin,
	# 	"radius": radius,
	# 	"theta_tw_1": theta_tw_1
	# }

	#geom_object = json.dumps(Geom, indent=2)
	#with open("geometry_file_test.json", "w") as geom_file:
	#	geom_file.write(geom_object)

	# Data to be written in parameters file
	# param = {
	# 	"Aircraft": "xyz",
	# 	"computational parameters": [
	# 		{
	# 			"d_psi": 1.0,
	# 			"iterations": 5,
	# 			"requested elements": 48,
	# 			"shed_history_angle": 45.0,
	# 			"wake_dist": [ 2.5 ],
	# 			"wopwop_rotor_indx" : 0
	# 		}
	# 	],
	# 	"flight conditions": [
	# 		{
	# 			"aoa": aoas,
	# 			"aoa_rotors": aoa_rotors,
	# 			"density": 1.125,
	# 			"omega": omega,
	# 			"sos": 343,
	# 			"V inf": V_inf
	# 		}
	# 	]
	# }

	density = base_param["flight conditions"][0]["density"]
	omegas = base_param["flight conditions"][0]["omegas"]
	sos = base_param["flight conditions"][0]["sos"]
	aoa_rotors = base_param["flight conditions"][0]["aoa_rotors"]
	base_param["flight conditions"].clear()

	for aoa in aoas:
		for V_inf in V_infs:
			for c_t in c_ts:
				flight_condition = {
					"aoa": aoa,
					"aoa_rotors": aoa_rotors,
					"density": density,
					"omegas": omegas,
					"sos": sos,
					"c_t": [c_t for i in range(len(omegas))],
					"V_inf": V_inf
				}

				base_param["flight conditions"].append(flight_condition)

	param_object = json.dumps(base_param, indent=2)
	with open("parameters_file_test.json", "w") as param_file:
		param_file.write(param_object)

	geom_object = json.dumps(aircrafts, indent=2)
	with open("geometry_test_file.json", "w") as geom_file:
		geom_file.write(geom_object)

if __name__ == "__main__":
	main()
