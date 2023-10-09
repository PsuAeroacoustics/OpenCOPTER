import sys
import os

sys.path.append(f'{os.path.dirname(os.path.realpath(__file__))}/dependencies/OpenCOPTER')
sys.path.append(f'{os.path.dirname(os.path.realpath(__file__))}/dependencies/wopwopd')

from libopencopter import *
from libwopwopd import *
import numpy as np
import math

RAD_TO_HZ = 0.1591549

# def generate_wopwop_namelist(R, origins, atmo, num_rotors, num_blades, omegas, dt, V_inf, iterations, aoa, output_path, actual_rotor_idx, t_min, t_max, nt):
# 	aircraft_cob = CB()
# 	aircraft_cob.Title = "Forward Velocity"
# 	aircraft_cob.TranslationType = TranslationType_known_function()
# 	aircraft_cob.AH = FVec3([0, 0, 0])
# 	aircraft_cob.VH = FVec3([V_inf, 0, 0])
# 	aircraft_cob.Y0 = FVec3([0, 0, 0])

# 	wopwop_aircraft = ContainerIn()
# 	wopwop_aircraft.Title = "HART"
# 	wopwop_aircraft.dTau = dt
# 	wopwop_aircraft.tauMax = dt*iterations
# 	wopwop_aircraft.cobs = [aircraft_cob]

# 	environment_constants = EnvironmentConstants()
# 	environment_constants.rho = atmo.density
# 	environment_constants.c = 343
# 	environment_constants.nu = atmo.dynamic_viscosity

# 	environment_in = EnvironmentIn()
# 	environment_in.pressureFolderName = "pressure"
# 	environment_in.SPLFolderName = "spl"
# 	environment_in.sigmaFolderName = "sigma"
# 	environment_in.debugLevel = 12
# 	environment_in.ASCIIOutputFlag = False
# 	environment_in.OASPLdBFlag = True
# 	environment_in.OASPLdBAFlag = False
# 	environment_in.spectrumFlag = True
# 	environment_in.SPLdBFlag = True
# 	environment_in.SPLdBAFlag = False
# 	environment_in.pressureGradient1AFlag = False
# 	environment_in.acousticPressureFlag = True
# 	environment_in.thicknessNoiseFlag = True
# 	environment_in.loadingNoiseFlag = True
# 	environment_in.totalNoiseFlag = True
# 	environment_in.sigmaFlag = False

# 	def build_blade_cntr(blade_idx):
# 		blade_cntr = ContainerIn()
# 		blade_cntr.Title = "blade "+str(blade_idx)+" containter"
# 		blade_cntr.patchGeometryFile = "blade_geom.dat"
# 		blade_cntr.patchLoadingFile = f"blade_{actual_rotor_idx}_{blade_idx}_loading.dat"

# 		blade_cob = CB()
# 		blade_cob.Title = "blade "+str(blade_idx)+" azimuth offset"
# 		blade_cob.AxisType = AxisType_time_independant()
# 		blade_cob.AngleType = AngleType_time_independant()
# 		blade_cob.AxisValue = FVec3([0, 0, 1])
# 		blade_cob.AngleValue = blade_idx*(2.0*math.pi/num_blades)
# 		blade_cob.Psi0 = blade_idx*(2.0*math.pi/num_blades)

# 		blade_cntr.cobs = [blade_cob]

# 		return blade_cntr

# 	def build_rotor_cntr(rotor_idx):
# 		rotor_cntr = ContainerIn()
# 		rotor_cntr.Title = "Rotor "+str(rotor_idx)

# 		rotor_offset_cb = CB()
# 		rotor_offset_cb.Title = "Rotor "+str(actual_rotor_idx)+" vertical offset"
# 		rotor_offset_cb.TranslationType = TranslationType_time_independant()
# 		rotor_offset_cb.TranslationValue = FVec3([0, 0, R[rotor_idx]*origins[rotor_idx][2]])

# 		rotor_aoa_cb = CB()
# 		rotor_aoa_cb.Title = "Rotor aoa"
# 		rotor_aoa_cb.AxisType = AxisType_time_independant()
# 		rotor_aoa_cb.AxisValue = FVec3([0, 1, 0])
# 		rotor_aoa_cb.AngleType = AngleType_time_independant()
# 		rotor_aoa_cb.AngleValue = -aoa

# 		rotor_rotation_cb = CB()
# 		rotor_rotation_cb.Title = "Rotor "+str(actual_rotor_idx)+" rotation"
# 		rotor_rotation_cb.AxisType = AxisType_time_independant()
# 		rotor_rotation_cb.AxisValue = FVec3([0, 0, 1])
# 		rotor_rotation_cb.AngleType = AngleType_known_function()
# 		rotor_rotation_cb.Omega = omegas[rotor_idx]
# 		rotor_rotation_cb.Rotation = True

# 		rotor_cntr.cobs = [rotor_offset_cb, rotor_aoa_cb, rotor_rotation_cb]

# 		rotor_cntr.children = [build_blade_cntr(blade_idx) for blade_idx in range(num_blades)]

# 		return rotor_cntr

# 	wopwop_aircraft.children = [build_rotor_cntr(rotor_idx) for rotor_idx in range(num_rotors)]
# 	RAD_TO_HZ = 1/(2*np.pi)
# 	blade_passing_freq = num_blades*abs(omegas[0])*RAD_TO_HZ
# 	lower_mid_freq = 6*blade_passing_freq
# 	upper_mid_freq = 25*blade_passing_freq

# 	#observer_nt = int(4*(1.0*(2.0*math.pi/abs(omega))/dt))
# 	#tMin = 0.372249
# 	observer = ObserverIn()
# 	observer.Title = "Mic array"
# 	observer.nt = nt# observer_nt
# 	observer.tMin = t_min #tMin
# 	observer.tMax = t_max # tMin + 2.0*(2.0*math.pi/abs(omega))

# 	#observer.xLoc = 16.577
# 	#observer.yLoc = 14.6972
# 	#observer.zLoc = -22.9814

# 	observer.nbTheta = 10
# 	observer.nbPsi = 10
# 	observer.thetaMin = 0
# 	observer.thetaMax = 2*math.pi
# 	observer.psiMin = -math.pi/2.0
# 	observer.psiMax = 0
# 	observer.radius = 20*R[0]
# 	observer.nbHarmonics = 30
# 	# observer.nbx = 23
# 	# observer.nby = 17
# 	# observer.nbz = 1
# 	# observer.xMin = -2.0*R[0]
# 	# observer.xMax = 2.0*R[0]
# 	# observer.yMin = -1.345*R[0]
# 	# observer.yMax = 1.345*R[0]
# 	# observer.zMin = -1.1075*R[0]
# 	# observer.zMax = -1.1075*R[0]

# 	observer.lowPassFrequency = upper_mid_freq
# 	observer.highPassFrequency = lower_mid_freq
# 	observer.cobs = [aircraft_cob]

# 	namelist = Namelist()
# 	namelist.environment_in = environment_in
# 	namelist.environment_constants = environment_constants
# 	namelist.observers = [observer]
# 	namelist.containers = [wopwop_aircraft]

# 	write_namelist(namelist, f"{output_path}/case.nam")

# 	return namelist

def generate_wopwop_namelist(R, origins, atmo, num_rotors, num_blades, omegas, dt, V_inf, iterations, aoa, actual_rotor_idx, t_min, t_max, nt, observer_config, acoustics_config, wopwop_data_path, sos, collective, full_system = False):
	aircraft_cob = CB()
	aircraft_cob.Title = "Forward Velocity"
	aircraft_cob.TranslationType = TranslationType_known_function()
	aircraft_cob.AH = FVec3([0, 0, 0])
	aircraft_cob.VH = FVec3([V_inf, 0, 0])
	aircraft_cob.Y0 = FVec3([0, 0, 0])

	aircraft_aoa_cb = CB()
	aircraft_aoa_cb.Title = "Aircraft aoa"
	aircraft_aoa_cb.AxisType = AxisType_time_independant()
	aircraft_aoa_cb.AxisValue = FVec3([0, 1, 0])
	aircraft_aoa_cb.AngleType = AngleType_time_independant()
	aircraft_aoa_cb.AngleValue = aoa

	wopwop_aircraft = ContainerIn()
	wopwop_aircraft.Title = "Aircraft"
	wopwop_aircraft.dTau = dt
	wopwop_aircraft.tauMax = dt*iterations
	wopwop_aircraft.cobs = [aircraft_cob, aircraft_aoa_cb]

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
	
	def build_blade_cntr(rotor_idx, blade_idx):
		blade_cntr = ContainerIn()
		blade_cntr.Title = "blade "+str(blade_idx)+" containter"
		blade_cntr.patchGeometryFile = f"../data/blade_geom_{rotor_idx}_{blade_idx}.dat"
		blade_cntr.patchLoadingFile = f"../data/blade_{rotor_idx if full_system else actual_rotor_idx}_{blade_idx}_loading.dat"

		blade_cob = CB()
		blade_cob.Title = "blade "+str(blade_idx)+" azimuth offset"
		blade_cob.AxisType = AxisType_time_independant()
		blade_cob.AngleType = AngleType_time_independant()
		blade_cob.AxisValue = FVec3([0, 0, 1])

		blade_collective_cob = CB()
		blade_collective_cob.Title = "blade "+str(blade_idx)+" pitch offset"
		blade_collective_cob.AxisType = AxisType_time_independant()
		blade_collective_cob.AngleType = AngleType_time_independant()
		blade_collective_cob.AxisValue = FVec3([1, 0, 0])

		if isinstance(num_blades, list):
			blade_cob.AngleValue = blade_idx*(2.0*math.pi/num_blades[rotor_idx])
			blade_cob.Psi0 = blade_idx*(2.0*math.pi/num_blades[rotor_idx])
		else:
			blade_cob.AngleValue = blade_idx*(2.0*math.pi/num_blades)
			blade_cob.Psi0 = blade_idx*(2.0*math.pi/num_blades)

		if omegas[rotor_idx] < 0:
			blade_collective_cob.AngleValue = math.pi - collective[rotor_idx]
		else:
			blade_collective_cob.AngleValue = collective[rotor_idx]


		blade_cntr.cobs = [blade_cob, blade_collective_cob]

		return blade_cntr

	def build_rotor_cntr(rotor_idx):
		rotor_cntr = ContainerIn()
		rotor_cntr.Title = "Rotor "+str(rotor_idx)

		rotor_offset_cb = CB()
		rotor_offset_cb.Title = "Rotor "+str(rotor_idx if full_system else actual_rotor_idx)+" vertical offset"
		rotor_offset_cb.TranslationType = TranslationType_time_independant()
		if full_system:
			rotor_offset_cb.TranslationValue = FVec3([R[rotor_idx]*origins[rotor_idx][0], R[rotor_idx]*origins[rotor_idx][1], R[rotor_idx]*origins[rotor_idx][2]])
		else:
			rotor_offset_cb.TranslationValue = FVec3([R[actual_rotor_idx]*origins[rotor_idx][0], R[actual_rotor_idx]*origins[rotor_idx][1], R[actual_rotor_idx]*origins[rotor_idx][2]])

		#rotor_aoa_cb = CB()
		#rotor_aoa_cb.Title = "Rotor aoa"
		#rotor_aoa_cb.AxisType = AxisType_time_independant()
		#rotor_aoa_cb.AxisValue = FVec3([0, 1, 0])
		#rotor_aoa_cb.AngleType = AngleType_time_independant()
		#rotor_aoa_cb.AngleValue = -aoa

		rotor_rotation_cb = CB()
		rotor_rotation_cb.Title = "Rotor "+str(rotor_idx if full_system else actual_rotor_idx)+" rotation"
		rotor_rotation_cb.AxisType = AxisType_time_independant()
		rotor_rotation_cb.AxisValue = FVec3([0, 0, 1])
		rotor_rotation_cb.AngleType = AngleType_known_function()
		rotor_rotation_cb.Omega = omegas[rotor_idx]
		rotor_rotation_cb.Rotation = True

		rotor_cntr.cobs = [rotor_offset_cb, rotor_rotation_cb]

		if isinstance(num_blades, list):
			rotor_cntr.children = [build_blade_cntr(rotor_idx, blade_idx) for blade_idx in range(num_blades[rotor_idx])]
		else:
			rotor_cntr.children = [build_blade_cntr(rotor_idx, blade_idx) for blade_idx in range(num_blades)]

		return rotor_cntr

	wopwop_aircraft.children = [build_rotor_cntr(rotor_idx) for rotor_idx in range(num_rotors)]

	if isinstance(num_blades, list):
		blade_passing_freq = min(num_blades)*abs(omegas[0])*RAD_TO_HZ
	else:
		blade_passing_freq = num_blades*abs(omegas[0])*RAD_TO_HZ
	
	observer = ObserverIn()

	observer.Title = "Mic array"
	observer.nt = nt
	observer.tMin = t_min
	observer.tMax = t_max

	if observer_config["type"] == "plane":
		if observer_config["radii_relative"]:
			observer.nbx = observer_config["nb"][0]*np.max(R)
			observer.nby = observer_config["nb"][1]*np.max(R)
			observer.nbz = observer_config["nb"][2]*np.max(R)
			observer.xMin = observer_config["min"][0]*np.max(R)
			observer.xMax = observer_config["max"][0]*np.max(R)
			observer.yMin = observer_config["min"][1]*np.max(R)
			observer.yMax = observer_config["max"][1]*np.max(R)
			observer.zMin = observer_config["min"][2]*np.max(R)
			observer.zMax = observer_config["max"][2]*np.max(R)
		else:
			observer.nbx = observer_config["nb"][0]
			observer.nby = observer_config["nb"][1]
			observer.nbz = observer_config["nb"][2]
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
			observer.radius = observer_config["radius"]*np.max(R)
		else:
			observer.radius = observer_config["radius"]

	elif observer_config["type"] == "points":
		with open(f"{wopwop_data_path}/observers.dat", 'w') as observer_file:
			observer_file.write(f"{''.join([f'{l} ' for l in observer_config['layout']])}\n")
			observer_file.write(f"{''.join([f'{l} ' for l in observer_config['x']])}\n")
			observer_file.write(f"{''.join([f'{l} ' for l in observer_config['y']])}\n")
			observer_file.write(f"{''.join([f'{l} ' for l in observer_config['z']])}\n")

			observer.fileName = f"../data/observers.dat"

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

	namelist = Namelist()
	namelist.environment_in = environment_in
	namelist.environment_constants = environment_constants
	namelist.observers = [observer]
	namelist.containers = [wopwop_aircraft]

	return namelist


def write_wopwop_geometry(airfoil_xsection, r, twist, R, real_chord, output_path, rotor_idx, blade_idx):
	print("Building wopwop geometry")

	#real_chord = (R/AR)*np.ones(len(r))
	blade_geom = generate_simple_constant_blade_geom(airfoil_xsection, r, twist, R, real_chord)
	lifting_line_geom = GeometryData(num_nodes = len(r))

	lifting_line_geom.set_x_nodes([R*_r for _r in r])
	lifting_line_geom.set_y_nodes(np.zeros(len(r), dtype=np.single))
	lifting_line_geom.set_z_nodes(np.zeros(len(r), dtype=np.single))

	lifting_line_geom.set_x_normals(np.zeros(len(r), dtype=np.single))
	lifting_line_geom.set_y_normals(np.zeros(len(r), dtype=np.single))
	lifting_line_geom.set_z_normals(np.ones(len(r), dtype=np.single))

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

	geom_file = create_geometry_file(wopwop_geom, f"{output_path}/blade_geom_{rotor_idx}_{blade_idx}.dat")
	append_geometry_data(geom_file, lifting_line_geom, 0)
	append_geometry_data(geom_file, blade_geom, 1)
	#append_geometry_data(geom_file, blade_geom, 0)
	
	close_geometry_file(geom_file)

def build_wopwop_loading(rotor_idx, blade_idx, iterations, r, airfoil_xsection, output_path):
    loading = AperiodicStructuredVectorLoadingFile(
		comment = "Blade loading",
		reference_frame = ReferenceFrame_patch_fixed(),
		data_alignment = DataAlignment_node_centered(),
		zone_headers = [
			AperiodicStructuredHeader(
				name = "Lift line loading",
				timesteps = iterations,
				i_max = len(r),
				j_max = 1,
				zone = 1,
				compute_thickness = False,
				has_data = True
			),
			AperiodicStructuredHeader(
				name = "Dummy blade loading",
				timesteps = iterations,
				i_max = len(r),
				j_max = len(airfoil_xsection),
				zone = 2,
				compute_thickness = True,
				has_data = False
			)
		]
	)
    
    loading_file = create_loading_file(loading, f"{output_path}/blade_{rotor_idx}_{blade_idx}_loading.dat")
    return loading_file