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

def generate_wopwop_namelist(R, origins, atmo, num_rotors, num_blades, omegas, dt, V_inf, iterations, aoa, data_path, actual_rotor_idx, t_min, t_max, nt, full_system = False):
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
	aircraft_aoa_cb.AngleValue = -aoa

	wopwop_aircraft = ContainerIn()
	wopwop_aircraft.Title = "HART"
	wopwop_aircraft.dTau = dt
	wopwop_aircraft.tauMax = dt*iterations
	wopwop_aircraft.cobs = [aircraft_cob, aircraft_aoa_cb]

	environment_constants = EnvironmentConstants()
	environment_constants.rho = atmo.density
	environment_constants.c = 343
	environment_constants.nu = atmo.dynamic_viscosity

	environment_in = EnvironmentIn()
	environment_in.pressureFolderName = "pressure"
	environment_in.SPLFolderName = "spl"
	environment_in.sigmaFolderName = "sigma"
	environment_in.debugLevel = 12
	environment_in.ASCIIOutputFlag = False
	environment_in.OASPLdBFlag = True
	environment_in.OASPLdBAFlag = False
	environment_in.spectrumFlag = True
	environment_in.SPLdBFlag = True
	environment_in.SPLdBAFlag = False
	environment_in.pressureGradient1AFlag = False
	environment_in.acousticPressureFlag = True
	environment_in.thicknessNoiseFlag = True
	environment_in.loadingNoiseFlag = True
	environment_in.totalNoiseFlag = True
	environment_in.sigmaFlag = False

	def build_blade_cntr(rotor_idx, blade_idx):
		blade_cntr = ContainerIn()
		blade_cntr.Title = "blade "+str(blade_idx)+" containter"
		blade_cntr.patchGeometryFile = f"../data/blade_geom.dat"
		blade_cntr.patchLoadingFile = f"../data/blade_{rotor_idx if full_system else actual_rotor_idx}_{blade_idx}_loading.dat"

		blade_cob = CB()
		blade_cob.Title = "blade "+str(blade_idx)+" azimuth offset"
		blade_cob.AxisType = AxisType_time_independant()
		blade_cob.AngleType = AngleType_time_independant()
		blade_cob.AxisValue = FVec3([0, 0, 1])
		blade_cob.AngleValue = blade_idx*(2.0*math.pi/num_blades)
		blade_cob.Psi0 = blade_idx*(2.0*math.pi/num_blades)

		blade_cntr.cobs = [blade_cob]

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

		rotor_cntr.children = [build_blade_cntr(rotor_idx, blade_idx) for blade_idx in range(num_blades)]

		return rotor_cntr

	wopwop_aircraft.children = [build_rotor_cntr(rotor_idx) for rotor_idx in range(num_rotors)]

	blade_passing_freq = num_blades*abs(omegas[0])*RAD_TO_HZ
	lower_mid_freq = 6*blade_passing_freq
	upper_mid_freq = 25*blade_passing_freq

	#observer_nt = int(4*(1.0*(2.0*math.pi/abs(omega))/dt))
	#tMin = 0.372249
	observer = ObserverIn()
	observer.Title = "Mic array"
	observer.nt = nt# observer_nt
	observer.tMin = t_min #tMin
	observer.tMax = t_max # tMin + 2.0*(2.0*math.pi/abs(omega))

	#observer.xLoc = 16.577
	#observer.yLoc = 14.6972
	#observer.zLoc = -22.9814

	observer.nbTheta = 20
	observer.nbPsi = 20
	observer.thetaMin = 0
	observer.thetaMax = 2*math.pi
	observer.psiMin = -math.pi/2.0
	observer.psiMax = 0
	if full_system:
		observer.radius = 20*np.max(R)
	else:
		observer.radius = 20*R[actual_rotor_idx]
	#observer.nbHarmonics = 30
	#observer.nbx = 17
	#observer.nby = 17
	#observer.nbz = 1
	#observer.xMin = -20.0*R
	#observer.xMax = 20.0*R
	#observer.yMin = -20*R
	#observer.yMax = 20*R
	#observer.zMin = -20*R
	#observer.zMax = -20*R

	observer.lowPassFrequency = upper_mid_freq
	observer.highPassFrequency = lower_mid_freq
	observer.cobs = [aircraft_cob]
	#observer.cobs = [aircraft_aoa_cb]

	namelist = Namelist()
	namelist.environment_in = environment_in
	namelist.environment_constants = environment_constants
	namelist.observers = [observer]
	namelist.containers = [wopwop_aircraft]

	#write_namelist(namelist, f"{output_path}/case.nam")

	return namelist


def write_wopwop_geometry(airfoil_xsection, r, twist, R, AR, output_path):
	print("Building wopwop geometry")

	real_chord = (R/AR)*np.ones(len(r))
	blade_geom = generate_simple_constant_blade_geom(airfoil_xsection, r, twist, 2, real_chord)
	lifting_line_geom = GeometryData(num_nodes = len(r))

	lifting_line_geom.set_x_nodes([R*_r for _r in r])
	lifting_line_geom.set_y_nodes(np.zeros(len(r), dtype=np.single))
	lifting_line_geom.set_z_nodes(np.zeros(len(r), dtype=np.single))

	lifting_line_geom.set_x_normals(np.zeros(len(r), dtype=np.single))
	lifting_line_geom.set_y_normals(np.zeros(len(r), dtype=np.single))
	lifting_line_geom.set_z_normals(np.ones(len(r), dtype=np.single))

	wopwop_geom = GeometryFile(
		comment = "Tilt rotor blade geometry",
		units = "Pa",
		data_alignment = DataAlignment_node_centered(), # node centered,
		zone_headers = [
			ConstantStructuredGeometryHeader(
				name = "blade",
				i_max = len(r),
				j_max = len(airfoil_xsection)
			),
			ConstantStructuredGeometryHeader(
				name = "lifting line",
				i_max = len(r),
				j_max = 1
			)
		]
	)

	geom_file = create_geometry_file(wopwop_geom, f"{output_path}/blade_geom.dat")
	append_geometry_data(geom_file, blade_geom, 0)
	append_geometry_data(geom_file, lifting_line_geom, 1)
	close_geometry_file(geom_file)

def build_wopwop_loading(rotor_idx, blade_idx, iterations, r, airfoil_xsection, output_path):
    loading = AperiodicStructuredVectorLoadingFile(
		comment = "HART rotor blade loading",
		reference_frame = ReferenceFrame_patch_fixed(),
		data_alignment = DataAlignment_node_centered(),
		zone_headers = [
			AperiodicStructuredHeader(
				name = "Dummy blade loading",
				timesteps = iterations,
				i_max = len(r),
				j_max = len(airfoil_xsection),
				zone = 1,
				compute_thickness = True,
				has_data = False
			),
			AperiodicStructuredHeader(
				name = "Lift line loading",
				timesteps = iterations,
				i_max = len(r),
				j_max = 1,
				zone = 2,
				compute_thickness = False,
				has_data = True
			),
		]
	)
    
    loading_file = create_loading_file(loading, f"{output_path}/blade_{rotor_idx}_{blade_idx}_loading.dat")
    return loading_file