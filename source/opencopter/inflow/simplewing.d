module opencopter.inflow.simplewing;

import opencopter.aircraft;
import opencopter.config;
import opencopter.inflow;
import opencopter.math;
import opencopter.math.blas;
import opencopter.math.lapacke;
import opencopter.memory;
import opencopter.wake;

import numd.calculus.integration.forwardeuler;
import numd.calculus.integration.rk4;
import numd.utility : linspace;

import std.array;
import std.algorithm;
import std.complex;
import std.conv;
import std.math;
import std.range;
import std.stdio;
import std.typecons;



class SimpleWingT(ArrayContainer AC = ArrayContainer.none) {

	VortexFilamentT!AC wing_filament;

	private Vec3 wing_origin;
	private double C_L;

	this(double _C_L, double V_inf, double V_tip, double y_0, double c, Vec3 _origin) {
		
		C_L = _C_L;
		wing_origin = _origin;
		wing_filament = VortexFilamentT!AC(chunk_size);

		debug writeln("wing_filament.chunks.length: ", wing_filament.chunks.length);
		foreach(c_idx; 0..chunk_size) {
			wing_filament.chunks[0].x[c_idx] = 0;//origin[0];
			wing_filament.chunks[0].y[c_idx] = 2*y_0*c_idx.to!double/((chunk_size - 1).to!double) - y_0;// + origin[1];
			wing_filament.chunks[0].z[c_idx] = 0;//origin[2];
			wing_filament.chunks[0].r_c[c_idx] = 0.5*c;
			wing_filament.chunks[0].gamma[c_idx] = V_inf*C_L/(2.0*V_tip);
		}

		debug writeln("C_L: ", C_L);
		debug writeln("wing_origin: ", wing_origin);
		debug writeln("V_inf: ", V_inf);
		debug writeln("V_tip: ", V_tip);
		debug writeln("c: ", c);
		debug writeln("y_0: ", y_0);
	}

	double get_average_inflow() {
		return 0;
	}

	void update(double C_T, ref RotorInputStateT!AC rotor, ref RotorStateT!AC rotor_state, double advance_ratio, double axial_advance_ratio, double dt) {
		update_impl(C_T, rotor, rotor_state, advance_ratio, axial_advance_ratio, dt);
	}

	void update(double C_T, RotorInputStateT!AC* rotor, RotorStateT!AC* rotor_state, double advance_ratio, double axial_advance_ratio, double dt) {
		update_impl(C_T, rotor, rotor_state, advance_ratio, axial_advance_ratio, dt);
	}

	private void update_impl(RIS, RS)(double C_T, auto ref RIS rotor, auto ref RS rotor_state, double advance_ratio, double axial_advance_ratio, double dt)
		if(is_rotor_input_state!RIS)
	{

	}

	double wake_skew() {
		return 0;
	}
	
	Chunk inflow_at(immutable Chunk x, immutable Chunk y, immutable Chunk z, immutable Chunk x_e, double angle_of_attack) {
		auto ind_vel = compute_filament_induced_velocities(wing_filament.chunks, x, y, z, 0);

		immutable Chunk neg_z = -2.0*PI*ind_vel.v_z[];
		debug writeln("neg_z: ", neg_z);
		return neg_z;
	}

	Vec3 origin() { return wing_origin;}
}
