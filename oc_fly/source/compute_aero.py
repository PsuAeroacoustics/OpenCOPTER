
import sys
import os

sys.path.append(f'{os.path.dirname(os.path.realpath(__file__))}../../dependencies/OpenCOPTER')
sys.path.append(f'{os.path.dirname(os.path.realpath(__file__))}../../dependencies/wopwopd')

from copy import deepcopy
from libopencopter import *
from libwopwopd import *
import wopwop_input_files_generator
import simulate_aircraft
from simulated_vehicle import SimulatedVehicle
import numpy as np
import math
import scipy.io
from scipy.interpolate import interp1d
from scipy.misc import derivative
import numpy as np

from os import path, makedirs

def fourier_motion(a: list[float], b: list[float], dt: float, frame: Frame, motion_axis: Vec3, azimuth_offset: float, azimuth: float):
	h = 0
	h_star = 0

	for idx in range(len(a)):
		cos = math.cos(float(idx)*(azimuth + azimuth_offset))
		sin = math.sin(float(idx)*(azimuth + azimuth_offset))

		h = h + a[idx]*cos + b[idx - 1]*sin

		if idx > 0:
			h_star = h_star + (-float(idx)*a[idx]*sin) + float(idx)*b[idx - 1]*cos

	frame.set_rotation(motion_axis, h)

def constant_motion(omega: float, dt: float, frame: Frame, motion_axis: Vec3, azimuth: float):
	delta_azimuth = omega*dt

	frame.rotate(motion_axis, delta_azimuth)

def build_blade(blade_object, requested_elements, geom_directory, R, frame):

	theta_tw_1 = blade_object['theta_tw']
	airfoil_descs = blade_object['airfoils']
	r_c = 0
	if "r_c" in blade_object:
		r_c = blade_object["r_c"]

	r = generate_radius_points(requested_elements, r_c)

	#print(f'r = {r}')
	elements = len(r)

	non_dim_length = 1.0 - r_c
	
	airfoils = []
	extents = []
	for airfoil_desc in airfoil_descs:
		type = airfoil_desc['type']
		extent = [0, 0]
		af_extent = airfoil_desc['extent']
		airfoil = None
		if type == 'aerodas':
			airfoil = create_aerodas_from_xfoil_polar(f"{geom_directory}/{airfoil_desc['xfoil_polar']}", airfoil_desc['thickness'])
		elif type == 'thinaf':
			airfoil = ThinAirfoil(0.0401)
		elif type == "C81":
			airfoil = load_c81_file(airfoil_desc['filename'])
		else:
			print(f"Unsupported airfoil type: {type}. Defaulting to ThinAirfoil theory")
			airfoil = ThinAirfoil(0.0401)

		extent[0] = int(round((1.0 - (1.0/math.pi)*math.acos((2.0*(af_extent[0] - r_c)/non_dim_length) - 1.0))*elements))
		extent[1] = int(round((1.0 - (1.0/math.pi)*math.acos((2.0*(af_extent[1] - r_c)/non_dim_length) - 1.0))*elements)) - 1

		extents.append(extent)
		airfoils.append(airfoil)

	blade_airfoil = BladeAirfoil(airfoils, extents)

	x = np.zeros(elements)

	linear_x = None
	linear_r = None
	linear_c = None
	if 'x' in blade_object:
		linear_x = blade_object['x']
		linear_r = np.linspace(r_c, 1.0, len(linear_x))
		f_x = interp1d(linear_r, linear_x)
		x = f_x(r)

	c = None
	if "AR" in blade_object:
		AR = blade_object["AR"]
		#c = ((1.0 - r_c)/AR)*np.ones(elements)
		c = (1.0/AR)*np.ones(elements)
		if linear_r is not None:
			#linear_c = ((1.0 - r_c)/AR)*np.ones(len(linear_r))
			linear_c = (1.0/AR)*np.ones(len(linear_r))
		else:
			linear_r = np.linspace(r_c, 1.0, elements)
			linear_x = np.zeros(len(linear_r))
			#linear_c = ((1.0 - r_c)/AR)*np.ones(len(linear_r))
			linear_c = (1.0/AR)*np.ones(len(linear_r))

	else:
		linear_c = blade_object['c']
		if linear_r is None:
			linear_x = np.zeros(len(linear_c))
			linear_r = np.linspace(r_c, 1.0, len(linear_c))

		f_c = interp1d(linear_r, linear_c)
		c = f_c(r)

	x_over_c = np.asarray(linear_x)/np.asarray(linear_c)
	
	f_x_over_c = interp1d(linear_r, x_over_c, bounds_error=False, fill_value='extrapolate')

	linear_x_over_c_p = [derivative(f_x_over_c, _r, 1.0e-12) for _r in linear_r[3:-3]]

	f_x_over_c_p = interp1d(linear_r[3:-3], linear_x_over_c_p, bounds_error=False, fill_value='extrapolate')

	xp = f_x_over_c_p(r)

	twist = np.asarray([(_r - 0.75)*theta_tw_1*(1/(1 - r_c))*(math.pi/180.0) for _r in generate_radius_points(requested_elements, r_c)])

	# Build the geom of the blades
	blade = BladeGeometry(
		num_elements = elements,
		azimuth_offset = 0,
		average_chord = 0,
		airfoil = blade_airfoil,
		r_c = r_c
	)
	
	sweep = sweep_from_quarter_chord(r, x)

	#print(f'twist: {[t*(180/math.pi) for t in twist]}')

	# Convert from linear array format to OpenCOPTER chunk format
	set_r(blade, r)
	set_xi(blade, x)
	set_xi_p(blade, xp)
	set_twist(blade, twist)
	set_chord(blade, c)
	set_sweep(blade, sweep)

	compute_blade_vectors(blade)

	blade.average_chord = R*np.sum(c)/len(c)
	blade.blade_length = R*(1.0 - r_c)
	blade.frame = frame

	# if omega < 0:
	#     # neg_twist = -twist
	#     naca0012_xsection = [[(xsection[0] - 0.5), -xsection[1]] for xsection in naca0012_xsection]
	#     wopwop_input_files_generator.write_wopwop_geometry(naca0012_xsection, r, twist, R, real_chord, wopwop_data_path, r_idx, b_idx)
	# else:
	#     naca0012_xsection = [[-(xsection[0] - 0.5), -xsection[1]] for xsection in naca0012_xsection]
	#     wopwop_input_files_generator.write_wopwop_geometry(naca0012_xsection, r, twist, R, real_chord, wopwop_data_path, r_idx, b_idx)
	
	return blade

def build_component(component_json, parent_frame, components_ref_dict, rotor_ref_dict, blade_ref_dict, components_dict, current_rotor_radius, requested_elements, geom_directory, motion_axis_dict, motion, trim_frame_name, trim_axis_dict, ref_count, did_deref):

	rotors = []
	blades = []

	# Dereference component if needed
	if "ref" in component_json:
		did_deref = True
		referenced_component_name = component_json["ref"]

		if referenced_component_name in components_ref_dict:
			component_json = components_ref_dict[referenced_component_name]["obj"]
			ref_count = components_ref_dict[referenced_component_name]["ref_count"]
			components_ref_dict[referenced_component_name]["ref_count"] = ref_count + 1
		elif referenced_component_name in rotor_ref_dict:
			component_json = rotor_ref_dict[referenced_component_name]["obj"]
			ref_count = rotor_ref_dict[referenced_component_name]["ref_count"]
			rotor_ref_dict[referenced_component_name]["ref_count"] = ref_count + 1
		elif referenced_component_name in blade_ref_dict:
			component_json = blade_ref_dict[referenced_component_name]["obj"]
			ref_count = blade_ref_dict[referenced_component_name]["ref_count"]
			blade_ref_dict[referenced_component_name]["ref_count"] = ref_count + 1

	name = component_json["name"]
	frame_type = component_json["type"]
	origin = Vec3(component_json["origin"]) if "origin" in component_json else Vec3([0.0, 0.0, 0.0])
	angle_axis = Vec3(component_json["axis"]) if "axis" in component_json else Vec3([1.0, 0.0, 0.0])
	angle = component_json["axis_angle"]*(math.pi/180.0) if "axis_angle" in component_json else 0.0

	actual_name = name
	if did_deref:
		actual_name = name + " " + str(ref_count)

	components_dict[actual_name] = {
		"name": actual_name,
		"frame_type": frame_type,
		"origin": origin,
		"angle_axis": angle_axis,
		"angle": angle
	}

	if trim_frame_name is not None and trim_frame_name == name:
		trim_axis_dict[actual_name] = angle_axis


	if motion is not None:
		matched_motions = list(filter(lambda x: x["frame"] == name, motion))
		
		if len(matched_motions) > 0:
			#print(f'matched_motions: {matched_motions}')
			if "axis_angle_function" not in component_json:
				raise Exception(f"No motion function defined for frame {name}")

			motion_axis_dict[actual_name] = (angle_axis, component_json["axis_angle_function"])
			#print(f'motion_axis_dict: {motion_axis_dict}')

	if did_deref:
		name = component_json["name"] + " " + str(ref_count)

	component_frame = Frame(angle_axis, angle, origin, parent_frame, name, frame_type)

	if frame_type == FrameType_rotor():

		current_rotor_radius = component_json["radius"]

	elif frame_type == FrameType_blade():

		if current_rotor_radius == None:
			raise Exception("Found blade component without a parent rotor")

		blades.append(build_blade(component_json, requested_elements, geom_directory, current_rotor_radius, component_frame))

	child_components = []

	if "children" in component_json:
		for child_component in component_json["children"]:
			built_components, built_rotors, built_blades = build_component(child_component, component_frame, components_ref_dict, rotor_ref_dict, blade_ref_dict, components_dict, current_rotor_radius, requested_elements, geom_directory, motion_axis_dict, motion, trim_frame_name, trim_axis_dict, ref_count, did_deref)

			rotors = rotors + built_rotors
			blades = blades + built_blades

			child_components.append(built_components)

			if frame_type == FrameType_rotor():
				if len(built_blades) > 1:
					raise Exception("Rotor has multiple blades on single attachment")

				built_blades[0].azimuth_offset = child_component["axis_angle"]*(math.pi/180.0)

	elif frame_type != FrameType_blade():
		raise Exception("Component has no children and is not of type 'blade'")
	
	if frame_type == FrameType_rotor():

		rotating_rotor_frame = Frame(angle_axis, angle, origin, component_frame, name, frame_type)
		rotating_rotor_frame.children = child_components

		component_frame.name = component_frame.name + " fixed"

		component_frame.children = [rotating_rotor_frame]

		rotor = build_rotor(blades, rotating_rotor_frame, current_rotor_radius)

		rotors.append(rotor)

		blades = []

	else:
		component_frame.children = child_components


	return component_frame, rotors, blades

def build_rotor(blades, frame, radius):

	rotor = RotorGeometry(
		num_blades = len(blades),
		origin = Vec3([0.0, 0.0, 0.0]),
		radius = radius,
		solidity = len(blades)*blades[0].average_chord/(math.pi*radius)
	)

	rotor.blades = blades
	rotor.frame = frame

	return rotor

def count_rotors(component, rotor_ref_dict):
	num_rotors = 0

	if "type" in component:
		if "rotor" == component["type"]:
			num_rotors = 1

	if "ref" in component:
		referenced_component_name = component["ref"]
		if referenced_component_name in rotor_ref_dict:
			if "children" in rotor_ref_dict[referenced_component_name]:
				for child in rotor_ref_dict[referenced_component_name]["children"]:
					num_rotors = num_rotors + count_rotors(child) + 1
	else:
		if "children" in component:
			for child in component["children"]:
				num_rotors = num_rotors + count_rotors(child, rotor_ref_dict)

	return num_rotors

def build_aircraft(geometry, requested_elements, geom_directory, motion, trim_frame_name):
	print("building aircraft")

	components_ref_dict = {}
	blade_ref_dict = {}
	rotor_ref_dict = {}
	components_dict = {}
	if "blades" in geometry:
		for blade_obj in geometry['blades']:
			blade_ref_dict[blade_obj["name"]] = {"obj": blade_obj, "ref_count": 0}

	if "components" in geometry:
		for component_obj in geometry["components"]:
			components_ref_dict[component_obj["name"]] = {"obj": component_obj, "ref_count": 0}

	if "rotors" in geometry:
		for rotor_obj in geometry["rotors"]:
			rotor_ref_dict[rotor_obj["name"]] = {"obj": rotor_obj, "ref_count": 0}

	num_rotors = 0
	for child in geometry["children"]:
		num_rotors = count_rotors(child, rotor_ref_dict) + num_rotors

	aircraft = Aircraft(num_rotors)

	motion_dict = {}
	trim_axis_dict = {}

	root_frame, rotors, _ = build_component(geometry, aircraft.root_frame, components_ref_dict, rotor_ref_dict, blade_ref_dict, components_dict, None, requested_elements, geom_directory, motion_dict, motion, trim_frame_name, trim_axis_dict, 0, False)

	aircraft.rotors = rotors

	aircraft.root_frame.name = root_frame.name
	aircraft.root_frame.children = root_frame.children
	aircraft.root_frame.set_frame_type(FrameType_aircraft())

	aircraft.root_frame.update(Mat4_identity())

	return aircraft, motion_dict, trim_axis_dict, components_dict

def compute_aero(log_file, args, output_base, do_compute, case):

	flight_condition = case.condition
	computational_parameters = case.computational_parameters
	results = case.results
	observer = case.observer
	acoustics = case.acoustics
	wake_history_length = case.wake_lengths
	rotorcraft_system = case.aircraft
	motion_vec_dict = case.motion_vec_dict
	trim_vec_dict = case.trim_vec_dict
	component_dict = case.component_dict

	num_rotors = rotorcraft_system.rotors.length()

	omegas = []
	for m in flight_condition['motion']:
		if ('blade_element_func' in m) and (m['blade_element_func'] == "rotation"):
			omegas.append(m['omega'])
	omegas = np.asarray(omegas)

	d_psi = computational_parameters["d_psi"]*math.pi/180.0
	dt = d_psi/np.max(np.abs(omegas))

	motion_lambdas = [[] for r_idx in range(rotorcraft_system.rotors.length())]
	wopwop_motion = {}

	if "motion" in flight_condition:
		for r_idx, rotor in enumerate(rotorcraft_system.rotors):

			def build_motion_lambdas(parent_frame, azimuth_offset):
				sub_motion_lambdas = []
				for child in parent_frame.children:
					for child_motion in flight_condition["motion"]:
						child_motion = deepcopy(child_motion)
						if (child.name == child_motion["frame"]) or (child.name[0:-2] == child_motion["frame"]):

							if motion_vec_dict[child.name][1] == "fourier":
								motion_lambda = lambda a, cos=child_motion["cos"], sin=child_motion["sin"], dt=dt, frame=child, vec=motion_vec_dict[child.name][0], azimuth_offset=azimuth_offset: fourier_motion(cos, sin, dt, frame, vec, azimuth_offset, a)

								wopwop_motion[child.name] = {"type": "fourier", "A": child_motion["cos"], "B": child_motion["sin"], "vector": motion_vec_dict[child.name][0]}

							elif motion_vec_dict[child.name][1] == "constant":
								omega = motion["omega"]
								motion_lambda = lambda a, omega=omega, dt=dt, frame=child, vec=motion_vec_dict[rotor_frame.name][0]: constant_motion(omega, dt, frame, vec, a)
								wopwop_motion[child.name] = {"type": "constant", "omega": omega, "vector": motion_vec_dict[child.name][0]}

							sub_motion_lambdas.append(deepcopy(motion_lambda))

					sub_motion_lambdas = sub_motion_lambdas + build_motion_lambdas(child, azimuth_offset)

				return sub_motion_lambdas
			
			rotor_frame = rotor.frame

			for motion in flight_condition["motion"]:
				motion = deepcopy(motion)
				if rotor_frame.name == motion["frame"]:
					if motion_vec_dict[rotor_frame.name][1] == "fourier":

						#motion_lambda = lambda a: fourier_motion(motion["cos"], motion["sin"], dt, rotor.frame, motion_vec_dict[rotor_frame.name][0], 0, a)
						motion_lambda = lambda a, cos=motion["cos"], sin=motion["sin"], dt=dt, frame=rotor_frame, vec=motion_vec_dict[rotor_frame.name][0], azimuth_offset=azimuth_offset: fourier_motion(cos, sin, dt, frame, vec, azimuth_offset, a)
						wopwop_motion[rotor_frame.name] = {"type": "fourier", "A": motion["cos"], "B": motion["sin"], "vector": motion_vec_dict[rotor_frame.name][0]}

					elif motion_vec_dict[rotor_frame.name][1] == "constant":
						print(f"Adding constant motion lambda for rotor {rotor.frame.name}")
						omega = motion["omega"]
						motion_lambda = lambda a, omega=omega, dt=dt, frame=rotor.frame, vec=motion_vec_dict[rotor_frame.name][0]: constant_motion(omega, dt, frame, vec, a)
						wopwop_motion[rotor_frame.name] = {"type": "constant", "omega": omega, "vector": motion_vec_dict[rotor_frame.name][0]}

					motion_lambdas[r_idx].append(deepcopy(motion_lambda))

			for b_idx, attach_frame in enumerate(rotor_frame.children):
				motion_lambdas[r_idx] = motion_lambdas[r_idx] + build_motion_lambdas(attach_frame, rotor.blades[b_idx].azimuth_offset)

	trim_lambdas = []

	def build_trim_lambda(in_frame):
		lambdas = []
		for frame_name, axis in trim_vec_dict.items():
			if in_frame.name == frame_name:
				lambdas.append(lambda a: in_frame.set_rotation(axis, a))

		for child in in_frame.children:
			lambdas = lambdas + build_trim_lambda(child)

		return lambdas

	for r_idx, rotor in enumerate(rotorcraft_system.rotors):
		blade_lambdas = []
		for b_idx in range(rotor.blades.length()):
			blade_lambdas = blade_lambdas + build_trim_lambda(rotor.frame.children[b_idx])

		trim_lambdas.append(blade_lambdas)

	density = flight_condition["density"]
	dynamic_viscosity = 18.03e-6

	shed_history_angle = np.asarray(computational_parameters["shed_history_angle"])*math.pi/180.0
	shed_release_angle = computational_parameters["shed_release_angle"]*math.pi/180.0

	max_omega = np.max(np.abs(omegas))
	rotor_ratios = np.round(max_omega/np.abs(omegas))

	shed_history = np.round(shed_history_angle/(shed_release_angle)).astype(dtype=np.int64).tolist()
	release_ratio = np.round(rotor_ratios*shed_release_angle/d_psi).astype(dtype=np.int64).tolist()
	#rotor_ratios = rotor_ratios.astype(dtype=np.int64).tolist()
	
	print(f'shed_history: {shed_history}, release_ratio: {release_ratio}')
	requested_elements = computational_parameters["spanwise_elements"]

	R = np.asarray([rotor.radius for rotor in rotorcraft_system.rotors])

	mus = flight_condition["V_inf"]/(np.abs(omegas)*R)

	log_file.write(f'mus: {mus}\n')
	log_file.write(f'omegas: {omegas}\n')

	collectives = None
	if "collective" in flight_condition:
		collectives = np.asarray(flight_condition["collective"])*np.pi/180
	else:
		collectives = 3.0*np.ones(num_rotors)*np.pi/180

	atmo = Atmosphere(density = density, dynamic_viscosity = dynamic_viscosity, speed_of_sound = flight_condition["sos"])
	
	wake_history_revs = np.ceil(np.asarray(wake_history_length)*d_psi/(2.0*math.pi)).astype(int)
	log_file.write(f'wake_history_revs: {wake_history_revs}\n')
	
	wopwop_data_path = f'{output_base}/acoustics/data'

	if not path.isdir(wopwop_data_path):
		makedirs(wopwop_data_path, exist_ok=True)

	log_file.write(f'wake_history_length: {wake_history_length}\n') # this must be a list of integers

	num_blades = [rotor.blades.length() for rotor in rotorcraft_system.rotors]

	r = generate_radius_points(requested_elements)
	elements = len(r)
	log_file.write(f"requested_elements: {requested_elements}, actual elements:  {elements}\n")

	# AircraftState is the top level container for holding the current
	# aerodynamic state of the tandem_system. It breaks down into rotors
	# the the individual blades. There is a series of functions provided
	# to turn internal state data into a linear array.
	rotorcraft_state = None
	if do_compute:
		rotorcraft_state = AircraftState(num_rotors, num_blades, elements, rotorcraft_system)
		rotorcraft_state.freestream = Vec4([flight_condition["V_inf"], 0, 0, 0])

	log_file.write(f"Freestream vel: {flight_condition['V_inf']} m/s\n")

	# Create and setup the input state. This would be the
	# sort of input a dynamics simulator might feed into
	# the aero model.
	rotorcraft_input_state = AircraftInputState(num_rotors, num_blades)

	rotorcraft_system.root_frame.rotate(Vec3([0, 1, 0]), flight_condition["aoa"]*(math.pi/180.0))

	for r_idx in range(num_rotors):

		rotorcraft_input_state.rotor_inputs[r_idx].angular_velocity = omegas[r_idx]
		rotorcraft_input_state.rotor_inputs[r_idx].angular_accel = 0
		rotorcraft_input_state.rotor_inputs[r_idx].azimuth = 0 # azimuths[r_idx]

		for b_idx in range(num_blades[r_idx]):
			rotorcraft_input_state.rotor_inputs[r_idx].r_0[b_idx] = computational_parameters['r_0']*rotorcraft_system.rotors[r_idx].blades[b_idx].average_chord/rotorcraft_system.rotors[r_idx].radius
			rotorcraft_input_state.rotor_inputs[r_idx].blade_flapping[b_idx] = 0
			rotorcraft_input_state.rotor_inputs[r_idx].blade_flapping_rate[b_idx] = 0
			rotorcraft_input_state.rotor_inputs[r_idx].blade_pitches[b_idx] = collectives[r_idx]

	print(f'num_blades: {num_blades}')
	
	#rotorcraft_inflows = [HuangPeters(4, 2, rotorcraft_system.rotors[r_idx], dt) for r_idx in range(num_rotors)]
	rotorcraft_inflows = [HuangPeters(4, 2, rotorcraft_system.rotors[r_idx], dt) if num_blades[r_idx] != 2 else HuangPeters(2, 2, rotorcraft_system.rotors[r_idx], dt) for r_idx in range(num_rotors)]
	#rotorcraft_inflows = [HuangPeters(5, 3, rotorcraft_system.rotors[r_idx], dt) if num_blades[r_idx] != 2 else HuangPeters(2, 2, rotorcraft_system.rotors[r_idx], dt) for r_idx in range(num_rotors)]
	#rotorcraft_inflows = [HuangPeters(6, 2, rotorcraft_system.rotors[r_idx], dt) if num_blades[r_idx] != 2 else HuangPeters(2, 2, rotorcraft_system.rotors[r_idx], dt) for r_idx in range(num_rotors)]

	a1 = 6.5e-5
	if "a1" in computational_parameters:
		a1 = computational_parameters['a1']

	# Setup the wake history. We need at minimum 2 timesteps worth of history for the update.
	# Increasing the history increases computation time with the current implementation
	log_file.write(f'wake_history_length: {wake_history_length}\n')
	rotor_wake_history = WakeHistory(num_rotors, num_blades, wake_history_length, 2, elements, shed_history, release_ratio, a1)
	
	vehicle = SimulatedVehicle(
		rotorcraft_system,
		rotorcraft_state,
		rotorcraft_input_state,
		rotorcraft_inflows,
		rotor_wake_history,
		motion_lambdas,
		trim_lambdas
	)

	rotorcraft_thrusts, rotorcraft_namelists, results_dictionary = simulate_aircraft.simulate_aircraft(
		log_file,
		vehicle,
		atmo,
		elements,
		args.ws,
		f'{output_base}/vtu',
		f'{output_base}/acoustics',
		do_compute,
		flight_condition,
		computational_parameters,
		observer,
		acoustics,
		wake_history_length,
		results,
		wopwop_motion
	)
	
	if do_compute:
		# results_dictionary = {}
		for r_idx in range(num_rotors):
			actual_wake_history = wake_history_length[r_idx] if wake_history_length[r_idx]%chunk_size() == 0 else wake_history_length[r_idx] + (chunk_size() - wake_history_length[r_idx]%chunk_size())
			wake_trajectories = np.zeros((num_blades[r_idx], 3, actual_wake_history))
			wake_core_sizes = np.zeros((num_blades[r_idx], actual_wake_history))

			for b_idx in range(num_blades[r_idx]):
				wake_trajectories[b_idx, 0, :] = get_wake_x_component(rotor_wake_history.history[0].rotor_wakes[r_idx].tip_vortices[b_idx])
				wake_trajectories[b_idx, 1, :] = get_wake_y_component(rotor_wake_history.history[0].rotor_wakes[r_idx].tip_vortices[b_idx])
				wake_trajectories[b_idx, 2, :] = get_wake_z_component(rotor_wake_history.history[0].rotor_wakes[r_idx].tip_vortices[b_idx])
				wake_core_sizes[b_idx,  :] = get_wake_r_c_component(rotor_wake_history.history[0].rotor_wakes[r_idx].tip_vortices[b_idx])

			results_dictionary[f'wake_{r_idx}_trajectory'] = wake_trajectories
			results_dictionary[f'wake_{r_idx}_core_size'] = wake_core_sizes
			
		results_dictionary["rotor_c_t"] = rotorcraft_thrusts
		results_dictionary["rotor_collectives"] = [rotorcraft_input_state.rotor_inputs[r_idx].blade_pitches[0] for r_idx in range(num_rotors)]
		results_dictionary["rotor_chis"] = [rotorcraft_inflows[r_idx].wake_skew() for r_idx in range(num_rotors)]

		scipy.io.savemat(f"{output_base}/results.mat", results_dictionary)

		if args.vtu_results:
			if results is not None:
				if 'inflow_slices' in results:
					for slice_idx, inflow_slice in enumerate(results["inflow_slices"]):
						res_x = inflow_slice["resolution"][0]
						res_y = inflow_slice["resolution"][1]
						res_z = inflow_slice["resolution"][2]

						deltas = Vec3([inflow_slice["slice_size"][0]/res_x, inflow_slice["slice_size"][1]/res_y, inflow_slice["slice_size"][2]/res_z])
						start = Vec3(inflow_slice["slice_start"])

						#aoa =  vehicle.input_state.rotor_inputs[0].angle_of_attack
						#write_inflow_vtu(f"{vtu_output_path}/../inflow_model_slice_    .vtu", vehicle.inflows, deltas, start, res_x, res_y, res_z, 0, omegas[0], vehicle.aircraft.rotors)
						write_inflow_vtu(f"{output_base}/inflow_model_slice_{slice_idx}.vtu", rotorcraft_inflows, deltas, start, res_x, res_y, res_z, 0, omegas, rotorcraft_system.rotors)

				if 'wake_slices' in results:
					for slice_idx, inflow_slice in enumerate(results["wake_slices"]):
						res_x = inflow_slice["resolution"][0]
						res_y = inflow_slice["resolution"][1]
						res_z = inflow_slice["resolution"][2]

						deltas = Vec3([inflow_slice["slice_size"][0]/res_x, inflow_slice["slice_size"][1]/res_y, inflow_slice["slice_size"][2]/res_z])
						start = Vec3(inflow_slice["slice_start"])

						write_wake_field_vtu(f"{output_base}/wake_field_slice_{slice_idx}.vtu", rotorcraft_state, rotor_wake_history.history[0], deltas, start, res_x, res_y, res_z)

	cases = []
	if (acoustics is not None) and (observer is not None):
		if (num_rotors > 1) and (args.fs is False):
			for r_idx, namelist in enumerate(rotorcraft_namelists[0:num_rotors]):
				wopwop_case_path = f'{output_base}/acoustics/rotor_{r_idx}/'

				if not path.isdir(wopwop_case_path):
					makedirs(wopwop_case_path, exist_ok=True)

				single_rotor_case = Casename()
				single_rotor_case.caseNameFile = "case.nam"
				single_rotor_case.globalFolderName = wopwop_case_path
				single_rotor_case.namelist = namelist

				cases.append(single_rotor_case)


		wopwop_case_path = f'{output_base}/acoustics/full_system/'

		if not path.isdir(wopwop_case_path):
			makedirs(wopwop_case_path, exist_ok=True)

		full_system_case = Casename()
		full_system_case.caseNameFile = "case.nam"
		full_system_case.globalFolderName = wopwop_case_path
		full_system_case.namelist = rotorcraft_namelists[-1]

		cases.append(full_system_case)

	return cases
