import sys
import os

sys.path.append(f'{os.path.dirname(os.path.realpath(__file__))}/dependencies/OpenCOPTER')
sys.path.append(f'{os.path.dirname(os.path.realpath(__file__))}/dependencies/wopwopd')

from libopencopter import *
from libwopwopd import *
import numpy as np
import math

from os import path, makedirs

RAD_TO_HZ = 0.1591549

# lol
def flatten_children(frame, termination_type):
	def _flatten_children(frame, termination_type):

		if (frame.get_frame_type() != termination_type):

			components = _flatten_children(frame.parent, termination_type)

			if (frame.parent.get_frame_type() == FrameType_rotor()):
				components.append(frame)
			else:
				components.append(frame.parent)
				
			return components

		if (frame.get_frame_type() == FrameType_rotor()) or (frame.get_frame_type() == FrameType_aircraft()):
			return []
		else:
			return frame
		
	flat_children = _flatten_children(frame, termination_type)
	flat_children = flat_children[1:len(flat_children)]
	return flat_children

def make_cb(frame, wopwop_motion):
	new_cb = CB()
	new_cb.Title = frame.name
	new_cb.TranslationType = TranslationType_time_independant()

	frame_pos = frame.local_position()
	new_cb.TranslationValue = FVec3([frame_pos[0], frame_pos[1], frame_pos[2]])

	angle = frame.local_rotation_angle()
	axis = frame.local_rotation_axis()

	new_cb.AxisType = AxisType_time_independant()

	if not math.isnan(angle):
		new_cb.AxisValue = FVec3([axis[0], axis[1], axis[2]])

	if frame.name in wopwop_motion:
		motion = wopwop_motion[frame.name]
		motion_type = motion["type"]

		motion_axis = wopwop_motion[frame.name]['vector']
		new_cb.AxisValue = FVec3([motion_axis[0], motion_axis[1], motion_axis[2]])

		if motion_type == "constant":
			new_cb.AngleType = AngleType_known_function()
			new_cb.Omega = motion["omega"]
			new_cb.Rotation = True
		if motion_type == "fourier":
			new_cb.AngleType = AngleType_periodic()
			A = motion["A"]
			B = motion["B"]
			new_cb.A0 = A[0]
			new_cb.A = [-a for a in A[1:len(A)]]
			new_cb.B = [-b for b in B]

	else:
		new_cb.AngleType = AngleType_time_independant()
		if not math.isnan(angle):
			new_cb.AngleValue = angle
			if frame.parent and frame.parent.get_frame_type() == FrameType_rotor():
				new_cb.Psi0 = (angle - math.pi)
	
	return new_cb

def build_blade_cntr(rotor, blade, environment_in, wopwop_motion):
	blade_cntr = ContainerIn()

	blade_cntr.Title = blade.frame.name + " container"

	rotor_name = rotor.frame.name.replace(' ', '_').replace('\t', '_').replace('\n', '_')
	blade_name = blade.frame.name.replace(' ', '_').replace('\t', '_').replace('\n', '_')

	#if environment_in.thicknessNoiseFlag:
	blade_cntr.patchGeometryFile = f"../data/{rotor_name}_{blade_name}_geometry.dat"
	
	if environment_in.loadingNoiseFlag:
		blade_cntr.patchLoadingFile = f"../data/{rotor_name}_{blade_name}_loading.dat"

	flat_frame_list = flatten_children(blade.frame, FrameType_rotor())

	cob_list = [make_cb(frame, wopwop_motion) for frame in flat_frame_list]
	cob_list.append(make_cb(blade.frame, wopwop_motion))
	blade_cntr.cobs = cob_list

	return blade_cntr

def build_rotor_cntr(rotor, environment_in, wopwop_motion):
	rotor_cntr = ContainerIn()
	rotor_cntr.Title = rotor.frame.name + " container"

	flat_frame_list = flatten_children(rotor.frame, FrameType_aircraft())

	rotor_cntr.cobs = [make_cb(frame, wopwop_motion) for frame in flat_frame_list]

	rotor_cntr.children = [build_blade_cntr(rotor, blade, environment_in, wopwop_motion) for blade in rotor.blades]

	return rotor_cntr
	
def generate_wopwop_namelist(atmo, dt, V_inf, iterations, aoa, t_min, t_max, nt, observer_config, acoustics_config, wopwop_data_path, sos, aircraft, rotors, wopwop_motion, ac_input, wopwop_case_path):

	aircraft_cob = CB()
	aircraft_cob.Title = "Forward Velocity"
	aircraft_cob.TranslationType = TranslationType_known_function()
	aircraft_cob.AH = FVec3([0, 0, 0])
	aircraft_cob.VH = FVec3([-V_inf, 0, 0])
	aircraft_cob.Y0 = FVec3([0, 0, 0])

	aircraft_frame_change = CB()
	aircraft_frame_change.Title = "Frame Change"
	aircraft_frame_change.AxisType = AxisType_time_independant()
	aircraft_frame_change.AxisValue = FVec3([0, 0, 1])
	aircraft_frame_change.AngleType = AngleType_time_independant()
	aircraft_frame_change.AngleValue = math.pi

	aircraft_aoa_cb = CB()
	aircraft_aoa_cb.Title = "Aircraft aoa"
	aircraft_aoa_cb.AxisType = AxisType_time_independant()
	aircraft_aoa_cb.AxisValue = FVec3([0, 1, 0])
	aircraft_aoa_cb.AngleType = AngleType_time_independant()
	aircraft_aoa_cb.AngleValue = -aoa

	wopwop_aircraft = ContainerIn()
	wopwop_aircraft.Title = aircraft.root_frame.name + " container"
	wopwop_aircraft.dTau = dt
	wopwop_aircraft.tauMax = dt*iterations
	wopwop_aircraft.cobs = [aircraft_cob, aircraft_frame_change, aircraft_aoa_cb]

	environment_constants = EnvironmentConstants()
	environment_constants.rho = atmo.density
	environment_constants.c = sos
	environment_constants.nu = atmo.dynamic_viscosity

	environment_in = EnvironmentIn()

	environment_in.pressureFolderName = acoustics_config["pressure_folder_name"] if "pressure_folder_name" in acoustics_config else "pressure"
	environment_in.SPLFolderName = acoustics_config["spl_folder_name"] if "spl_folder_name" in acoustics_config else "spl"
	environment_in.sigmaFolderName = acoustics_config["sigma_folder_name"] if "sigma_folder_name" in acoustics_config else "sigma"
	environment_in.debugLevel = acoustics_config["debug_level"] if "debug_level" in acoustics_config else 12
	environment_in.ASCIIOutputFlag = acoustics_config["ascii_output_flag"] if "ascii_output_flag" in acoustics_config else False
	environment_in.OASPLdBFlag = acoustics_config["oaspl_db_flag"] if "oaspl_db_flag" in acoustics_config else True
	environment_in.OASPLdBAFlag = acoustics_config["oaspl_dba_flag"] if "oaspl_dba_flag" in acoustics_config else False
	environment_in.spectrumFlag = acoustics_config["spectrum_flag"] if "spectrum_flag" in acoustics_config else False
	environment_in.SPLdBFlag = acoustics_config["spl_db_flag"] if "spl_db_flag" in acoustics_config else True
	environment_in.SPLdBAFlag = acoustics_config["spl_dba_flag"] if "spl_dba_flag" in acoustics_config else False
	environment_in.pressureGradient1AFlag = acoustics_config["pressure_gradient_1a_flag"] if "pressure_gradient_1a_flag" in acoustics_config else False
	environment_in.acousticPressureFlag = acoustics_config["acoustic_pressure_flag"] if "acoustic_pressure_flag" in acoustics_config else True
	environment_in.thicknessNoiseFlag = acoustics_config["thickness_noise_flag"] if "thickness_noise_flag" in acoustics_config else True
	environment_in.loadingNoiseFlag = acoustics_config["loading_noise_flag"] if "loading_noise_flag" in acoustics_config else True
	environment_in.totalNoiseFlag = acoustics_config["total_noise_flag"] if "total_noise_flag" in acoustics_config else True
	environment_in.sigmaFlag = acoustics_config["sigma_flag"] if "sigma_flag" in acoustics_config else False
	environment_in.loadingNoiseSigmaFlag = acoustics_config["loading_noise_sigma_flag"] if "loading_noise_sigma_flag" in acoustics_config else False
	environment_in.thicknessNoiseSigmaFlag = acoustics_config["thickness_noise_sigma_flag"] if "thickness_noise_sigma_flag" in acoustics_config else False
	environment_in.totalNoiseSigmaFlag = acoustics_config["total_noise_sigma_flag"] if "total_noise_sigma_flag" in acoustics_config else False
	environment_in.normalSigmaFlag = acoustics_config["normal_sigma_flag"] if "normal_sigma_flag" in acoustics_config else False
	environment_in.machSigmaFlag = acoustics_config["mach_sigma_flag"] if "mach_sigma_flag" in acoustics_config else False
	environment_in.observerSigmaFlag = acoustics_config["observer_sigma_flag"] if "observer_sigma_flag" in acoustics_config else False
	environment_in.velocitySigmaFlag = acoustics_config["velocity_sigma_flag"] if "velocity_sigma_flag" in acoustics_config else False
	environment_in.accelerationSigmaFlag = acoustics_config["acceleration_sigma_flag"] if "acceleration_sigma_flag" in acoustics_config else False
	environment_in.densitySigmaFlag = acoustics_config["density_sigma_flag"] if "density_sigma_flag" in acoustics_config else False
	environment_in.momentumSigmaFlag = acoustics_config["momentum_sigma_flag"] if "momentum_sigma_flag" in acoustics_config else False
	environment_in.pressureSigmaFlag = acoustics_config["pressure_sigma_flag"] if "pressure_sigma_flag" in acoustics_config else False
	environment_in.loadingSigmaFlag = acoustics_config["loading_sigma_flag"] if "loading_sigma_flag" in acoustics_config else False
	environment_in.areaSigmaFlag = acoustics_config["area_sigma_flag"] if "area_sigma_flag" in acoustics_config else False
	environment_in.MdotrSigmaFlag = acoustics_config["mdotr_sigma_flag"] if "mdotr_sigma_flag" in acoustics_config else False
	environment_in.iblankSigmaFlag = acoustics_config["iblank_sigma_flag"] if "iblank_sigma_flag" in acoustics_config else False
	
	wopwop_aircraft.children = [build_rotor_cntr(rotor, environment_in, wopwop_motion) for rotor in rotors]

	R = 1
	num_blades = 1
	ref_omega = 1
	if ("reference_rotor" not in observer_config):
		raise Exception("Required key 'reference_rotor' not found in observer config file")
	elif ("radii_relative" in observer_config) and observer_config["radii_relative"]:
		R = aircraft.rotors[observer_config['reference_rotor']].radius
	
	num_blades = aircraft.rotors[observer_config['reference_rotor']].blades.length()
	ref_omega = ac_input.rotor_inputs[observer_config['reference_rotor']].angular_velocity

	blade_passing_freq = num_blades*abs(ref_omega)*RAD_TO_HZ
	
	observer = ObserverIn()

	observer.Title = "Mic array"
	observer.nt = nt
	observer.tMin = t_min
	observer.tMax = t_max

	if observer_config["type"] == "plane":
		observer.nbx = observer_config["nb"][0]
		observer.nby = observer_config["nb"][1]
		observer.nbz = observer_config["nb"][2]
		if observer_config["radii_relative"]:
			observer.xMin = observer_config["min"][0]*R
			observer.xMax = observer_config["max"][0]*R
			observer.yMin = observer_config["min"][1]*R
			observer.yMax = observer_config["max"][1]*R
			observer.zMin = observer_config["min"][2]*R
			observer.zMax = observer_config["max"][2]*R
		else:
			observer.xMin = observer_config["min"][0]
			observer.xMax = observer_config["max"][0]
			observer.yMin = observer_config["min"][1]
			observer.yMax = observer_config["max"][1]
			observer.zMin = observer_config["min"][2]
			observer.zMax = observer_config["max"][2]

	elif observer_config["type"] == "sphere":
		observer.nbTheta = observer_config["nb_theta"]
		observer.nbPsi = observer_config["nb_psi"]
		observer.thetaMin = observer_config["theta_min"]*(math.pi/180.0)
		observer.thetaMax = observer_config["theta_max"]*(math.pi/180.0)
		observer.psiMin = observer_config["psi_min"]*(math.pi/180.0)
		observer.psiMax = observer_config["psi_max"]*(math.pi/180.0)
		if observer_config["radii_relative"]:
			observer.radius = observer_config["radius"]*R
		else:
			observer.radius = observer_config["radius"]

	elif observer_config["type"] == "points":
		with open(f"{wopwop_data_path}/observers.dat", 'w') as observer_file:
			observer_file.write(f"{''.join([f'{l} ' for l in observer_config['layout']])}\n")
			observer_file.write(f"{''.join([f'{l} ' for l in observer_config['x']])}\n")
			observer_file.write(f"{''.join([f'{l} ' for l in observer_config['y']])}\n")
			observer_file.write(f"{''.join([f'{l} ' for l in observer_config['z']])}\n")

			observer.fileName = f"../data/observers.dat"

	elif observer_config["type"] == "single_point":
		if observer_config["radii_relative"]:
			observer.xLoc = observer_config["x"]*R
			observer.yLoc = observer_config["y"]*R
			observer.zLoc = observer_config["z"]*R
		else:
			observer.xLoc = observer_config["x"]
			observer.yLoc = observer_config["y"]
			observer.zLoc = observer_config["z"]

	if "low_pass_cutoff" in observer_config:
		if observer_config["cutoff_units"] == "bpf":
			observer.lowPassFrequency = observer_config["low_pass_cutoff"]*blade_passing_freq
		else:
			observer.lowPassFrequency = observer_config["low_pass_cutoff"]

	if "high_pass_cutoff" in observer_config:
		if observer_config["cutoff_units"] == "bpf":
			observer.highPassFrequency = observer_config["high_pass_cutoff"]*blade_passing_freq
		else:
			observer.highPassFrequency = observer_config["high_pass_cutoff"]

	observer.cobs = [aircraft_cob]

	if "frequency_ranges" in observer_config:

		if not path.isdir(f'{wopwop_case_path}/segmentProcess'):
			makedirs(f'{wopwop_case_path}/segmentProcess', exist_ok=True)

		range_list = []
		for frequency_range in observer_config["frequency_ranges"]:
			range_in = RangeIn()
			range_in.Title = frequency_range["name"]
			
			if "high_pass_cutoff" in frequency_range:
				if frequency_range["cutoff_units"] == "bpf":
					range_in.minFrequency = frequency_range["high_pass_cutoff"]*blade_passing_freq
				else:
					range_in.minFrequency = frequency_range["high_pass_cutoff"]

			if "low_pass_cutoff" in frequency_range:
				if frequency_range["cutoff_units"] == "bpf":
					range_in.maxFrequency = frequency_range["low_pass_cutoff"]*blade_passing_freq
				else:
					range_in.maxFrequency = frequency_range["low_pass_cutoff"]

			range_list.append(range_in)
		
		observer.ranges = range_list

	namelist = Namelist()
	namelist.environment_in = environment_in
	namelist.environment_constants = environment_constants
	namelist.observers = [observer]
	namelist.containers = [wopwop_aircraft]

	return namelist


def write_wopwop_geometry(airfoil_xsection, output_path, rotor, blade, include_thickness):
	print("Building wopwop geometry")

	R = rotor.radius
	r = get_r(blade)
	r = [_r - blade.r_c for _r in r]
	twist = get_twist(blade)
	real_chord = [c*R for c in get_chord(blade)]

	#real_chord = (R/AR)*np.ones(len(r))
	blade_geom = generate_simple_constant_blade_geom(airfoil_xsection, r, twist, R, real_chord)
	lifting_line_geom = GeometryData(num_nodes = len(r))

	lifting_line_geom.set_x_nodes([R*_r for _r in r])
	lifting_line_geom.set_y_nodes(np.zeros(len(r), dtype=np.single))
	lifting_line_geom.set_z_nodes(np.zeros(len(r), dtype=np.single))

	lifting_line_geom.set_x_normals(np.zeros(len(r), dtype=np.single))
	lifting_line_geom.set_y_normals(np.zeros(len(r), dtype=np.single))
	lifting_line_geom.set_z_normals(np.ones(len(r), dtype=np.single))

	rotor_name = rotor.frame.name.replace(' ', '_').replace('\t', '_').replace('\n', '_')
	blade_name = blade.frame.name.replace(' ', '_').replace('\t', '_').replace('\n', '_')

	if not include_thickness:
		wopwop_geom = GeometryFile(
			comment = "Blade geometry",
			units = "Pa",
			data_alignment = DataAlignment_node_centered(), # node centered,
			zone_headers = [
				ConstantStructuredGeometryHeader(
					name = "lifting line",
					i_max = len(r),
					j_max = 1
				),
				ConstantStructuredGeometryHeader(
					name = "blade",
					i_max = len(r),
					j_max = len(airfoil_xsection)
				)
			]
		)

		geom_file = create_geometry_file(wopwop_geom, f"{output_path}/{rotor_name}_{blade_name}_geometry.dat")
		
		append_geometry_data(geom_file, lifting_line_geom, 0)
		append_geometry_data(geom_file, blade_geom, 1)
	
	else:
		wopwop_geom = GeometryFile(
			comment = "Blade geometry",
			units = "Pa",
			data_alignment = DataAlignment_node_centered(), # node centered,
			zone_headers = [
				ConstantStructuredGeometryHeader(
					name = "lifting line",
					i_max = len(r),
					j_max = 1
				)
			]
		)

		geom_file = create_geometry_file(wopwop_geom, f"{output_path}/{rotor_name}_{blade_name}_geometry.dat")

		append_geometry_data(geom_file, lifting_line_geom, 0)


	close_geometry_file(geom_file)

def build_wopwop_loading(rotor, blade, iterations, airfoil_xsection, output_path):
	loading = AperiodicStructuredVectorLoadingFile(
		comment = "Blade loading",
		reference_frame = ReferenceFrame_patch_fixed(),
		data_alignment = DataAlignment_node_centered(),
		zone_headers = [
			AperiodicStructuredHeader(
				name = "Lift line loading",
				timesteps = iterations,
				i_max = len(get_r(blade)),
				j_max = 1,
				zone = 1,
				compute_thickness = False,
				has_data = True
			),
			AperiodicStructuredHeader(
				name = "Dummy blade loading",
				timesteps = iterations,
				i_max = len(get_r(blade)),
				j_max = len(airfoil_xsection),
				zone = 2,
				compute_thickness = True,
				has_data = False
			)
		]
	)
	
	rotor_name = rotor.frame.name.replace(' ', '_').replace('\t', '_').replace('\n', '_')
	blade_name = blade.frame.name.replace(' ', '_').replace('\t', '_').replace('\n', '_')
	loading_file = create_loading_file(loading, f"{output_path}/{rotor_name}_{blade_name}_loading.dat")
	return loading_file