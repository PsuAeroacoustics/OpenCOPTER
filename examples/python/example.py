#!/usr/bin/env python3
import sys
import os

# Tell python where it will find the libopencopter module
sys.path.append(f'{os.path.dirname(os.path.realpath(__file__))}/../../')

from libopencopter import *
import numpy as np
import math
import time

import matplotlib.pyplot as plt

if __name__ == "__main__":

	iterations = 5400
	#iterations = 900
	wake_history_length = 1*1024

	requested_elements = 45

	num_rotors = 1
	num_blades = 4
	R = 2

	theta_75 = 2.0*(math.pi/180.0)
	
	density = 1.125
	omega = 109.12
	AR = 16.5
	sos = 343

	V_inf = 32.9
	aoa = 6.5*(math.pi/180.0)
	theta_tw_1 = -6.5*(math.pi/180.0)

	d_psi = 1.0
	dt = d_psi*(math.pi/180.0)/math.fabs(omega)

	shed_history_angle = 45.0
	shed_history = round(shed_history_angle/d_psi)

	d_azimuth = 2.0*math.pi/num_blades

	# Root cutout
	r_c = 0.22

	# Create the spanwise distributions for our blades
	# We will later apply these to opencopters internal
	# data layout. r and c arrays are normalized by the
	# rotor radius
	r = generate_radius_points(requested_elements, r_c)
	elements = len(r)
	print("requested_elements: ", requested_elements, " actual elements: ", elements)
	c = (1.0/AR)*np.ones(elements)
	alpha_0 = np.zeros(elements)
	xi = np.zeros(elements)
	xi_p = np.zeros(elements)
	twist = [(_r - 0.75)*theta_tw_1 for _r in r]
	C_l_alpha = 2.0*math.pi*np.ones(elements)
	sweep = np.zeros(elements)

	atmo = Atmosphere(density = 1.125, dynamic_viscosity = 18.03e-6, speed_of_sound = 343)
	
	# Create our outer aircraft geometry container
	aircraft = Aircraft(num_rotors)
	
	# The origins of our 2 rotors
	origins = [Vec3([0, 0, 0]), Vec3([0, -2.5, 1])]

	def build_blade(b_idx, parent_frame):

		airfoil = create_aerodas_from_xfoil_polar(f"./NACA23012mod_1000000_polar.dat", 0.12)
		extent = [0, 47]
		blade_airfoil = BladeAirfoil([airfoil], [extent])

		# Build the geom of the blades
		blade = BladeGeometry(
			num_elements = elements,
			azimuth_offset = 0,
			average_chord = R*np.sum(c)/len(c),
			airfoil = blade_airfoil,
			r_c = r_c
		)

		blade.blade_length = R*(1.0 - r_c)

		# Convert from linear array format to opencopter chunk format
		set_r(blade, r)
		set_twist(blade, twist)
		set_alpha_0(blade, alpha_0)
		set_C_l_alpha(blade, C_l_alpha)
		set_chord(blade, c)
		set_sweep(blade, sweep)
		set_xi(blade, xi)
		set_xi_p(blade, xi_p)

		compute_blade_vectors(blade)

		azimuth_offset_frame = Frame(Vec3([0, 0, 1]), b_idx*d_azimuth, Vec3([0, 0, 0]), parent_frame, f'blade_{b_idx}_azimuth', FrameType_connection())
		root_cutout_frame = Frame(Vec3([1, 0, 0]), 0, Vec3([R*r_c, 0, 0]), azimuth_offset_frame, f'blade_{b_idx}_cutout', FrameType_connection())
		blade_frame = Frame(Vec3([1, 0, 0]), 0, Vec3([0, 0, 0]), root_cutout_frame, f'blade_{b_idx}', FrameType_blade())

		root_cutout_frame.children = [blade_frame]
		azimuth_offset_frame.children = [root_cutout_frame]

		blade.frame = blade_frame

		return blade

	
	def build_rotor(r_idx, parent_frame):
		rotor = RotorGeometry(
			num_blades = num_blades,
			origin = origins[r_idx],
			radius = R,
			solidity = 0 # we set this later because we don't have the average blade chord
		)

		rotor_frame = Frame(Vec3([1, 0, 0]), 0, Vec3([0, 0, 0]), parent_frame, f'rotor_{r_idx}', FrameType_rotor())

		rotor.blades = [build_blade(b_idx, rotor_frame) for b_idx in range(num_blades)]
		print(f'rotor.blades length: {rotor.blades.length()}')
		rotor.solidity = num_blades*rotor.blades[0].average_chord/(math.pi*rotor.radius)

		rotor_frame.children = [b.frame.parent.parent for b in rotor.blades]
		rotor.frame = rotor_frame

		return rotor

	aircraft.rotors = [build_rotor(r_idx, aircraft.root_frame) for r_idx in range(num_rotors)]
	
	aircraft.root_frame.name = "Aircraft frame"
	aircraft.root_frame.children = [r.frame for r in aircraft.rotors]
	aircraft.root_frame.set_frame_type(FrameType_aircraft())
	
	aircraft.root_frame.set_rotation(Vec3([0, 1, 0]), -aoa)

	aircraft.root_frame.update(Mat4_identity())

	# AircraftState is the top level container for holding the current
	# aerodynamic state of the aircraft. It breaks down into rotors
	# then the individual blades. There is a series of functions provided
	# to turn internal state data into a linear array.
	ac_state = AircraftState(num_rotors, [num_blades], elements, aircraft)

	print("Freestream vel: ", V_inf, " m/s")

	# Create and setup the input state. This would be the
	# sort of input a dynamics simulator might feed into
	# the aero model.
	ac_input_state = AircraftInputState(num_rotors, [num_blades])

	ac_state.freestream = Vec4([V_inf, 0, 0, 0])

	ac_input_state.rotor_inputs[0].angular_velocity = omega
	ac_input_state.rotor_inputs[0].angular_accel = 0
	ac_input_state.rotor_inputs[0].azimuth = 0
	for b_idx in range(num_blades):
		ac_input_state.rotor_inputs[0].r_0[b_idx] = 1.0*aircraft.rotors[0].blades[b_idx].average_chord
		ac_input_state.rotor_inputs[0].blade_flapping[b_idx] = 0
		ac_input_state.rotor_inputs[0].blade_flapping_rate[b_idx] = 0
		ac_input_state.rotor_inputs[0].blade_pitches[b_idx] = theta_75

	if num_rotors == 2:
		ac_input_state.rotor_inputs[1].angular_velocity = -omega # Negative indicates it rotates the other direction
		ac_input_state.rotor_inputs[1].angular_accel = 0
		ac_input_state.rotor_inputs[1].azimuth = 0

		for b_idx in range(num_blades):
			ac_input_state.rotor_inputs[1].r_0[b_idx] = 1.0*aircraft.rotors[1].blades[b_idx].average_chord
			ac_input_state.rotor_inputs[1].blade_flapping[b_idx] = 0
			ac_input_state.rotor_inputs[1].blade_flapping_rate[b_idx] = 0
			ac_input_state.rotor_inputs[1].blade_pitches[b_idx] = theta_75

	r_0 = ac_input_state.rotor_inputs[0].r_0[b_idx]
	print("r_0: ", ac_input_state.rotor_inputs[0].r_0)
	inflows = [HuangPeters(4, 2, aircraft.rotors[r_idx], dt) for r_idx in range(num_rotors)]

	# Setup the wake history. We need at minimum 2 timesteps worth of history for the update.
	# Increasing the history increases computation time with the current implementation
	wake_history = WakeHistory(num_rotors, [num_blades], [wake_history_length], 2, elements, [shed_history], [1])

	vtk_rotor = build_base_vtu_rotor(aircraft.rotors[0])
	vtk_wake = build_base_vtu_wake(wake_history.history[0])

	start_time = time.monotonic()

	for iteration in range(iterations):
		if iteration % 100 == 0:
			now = time.monotonic()
			print(now - start_time, ": iteration: ", iteration)
			start_time = now

		for r_idx in range(ac_input_state.rotor_inputs.length()):
			ac_input_state.rotor_inputs[r_idx].azimuth += ac_input_state.rotor_inputs[r_idx].angular_velocity*dt + ac_input_state.rotor_inputs[r_idx].angular_accel*dt*dt

			if ac_input_state.rotor_inputs[r_idx].azimuth > 2.0*math.pi:
				ac_input_state.rotor_inputs[r_idx].azimuth = math.fmod(ac_input_state.rotor_inputs[r_idx].azimuth, 2.0*math.pi)

			aircraft.rotors[r_idx].frame.set_rotation(Vec3([0, 0, 1]), ac_input_state.rotor_inputs[r_idx].azimuth)

		step(ac_state, aircraft, ac_input_state, inflows, wake_history, atmo, iteration, dt)

		if iteration > (iterations - 360):
			write_rotor_vtu("rotor", iteration, 0, vtk_rotor, ac_state.rotor_states[0], ac_input_state.rotor_inputs[0], aircraft.rotors[0])
			write_wake_vtu("wake", iteration, vtk_wake, wake_history.history[0])

	print("rotor 0 C_T: ", ac_state.rotor_states[0].C_T)

	max_dim_1 = -math.inf
	min_dim_1 = math.inf

	max_dim_2 = -math.inf
	min_dim_2 = math.inf

	# Iterating directly on the collections like this is read only
	for r_idx, rotor in enumerate(ac_state.rotor_states):
		for b_idx, blade in enumerate(rotor.blade_states):
			x = get_wake_x_component(wake_history.history[0].rotor_wakes[r_idx].tip_vortices[b_idx])
			z = get_wake_z_component(wake_history.history[0].rotor_wakes[r_idx].tip_vortices[b_idx])

			max_dim_2 = max(max_dim_2, max(z))
			min_dim_2 = min(min_dim_2, min(z))

			max_dim_1 = max(max_dim_1, max(x))
			min_dim_1 = min(min_dim_1, min(x))
			plt.plot(x, z, linewidth=0.5)

			x_r = math.cos(blade.azimuth)
			x_b = aircraft.rotors[r_idx].origin[0] + x_r*math.cos(aoa)
			z_b = aircraft.rotors[r_idx].origin[2] - x_r*math.sin(aoa)

			plt.plot([aircraft.rotors[r_idx].origin[0], x_b], [aircraft.rotors[r_idx].origin[2], z_b], "k-", linewidth=1.5)

	plt.axis("square")

	span_1 = max_dim_1 - min_dim_1
	span_2 = max_dim_2 - min_dim_2
	plt.xlim(left=min_dim_1 - .1*span_1, right= max_dim_1 + .1*span_1)
	plt.ylim(bottom=min_dim_2 - .9*span_2, top=max_dim_2 + .9*span_2)
	plt.show()
