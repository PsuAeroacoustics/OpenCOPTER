module opencopter.inflow;

import opencopter.aircraft;
import opencopter.memory;
import opencopter.math;
import opencopter.wake;

public import opencopter.inflow.beddos;
public import opencopter.inflow.huangpeters;
public import opencopter.inflow.simplewing;
public import opencopter.inflow.winginflow;

import numd.linearalgebra.matrix;

import std.stdio;
import std.range;
import std.algorithm;
import std.math;
import std.conv;
enum Direction {
	clockwise,
	counter_clockwise
}

alias Inflow = InflowT!(ArrayContainer.none);


class InflowT(ArrayContainer AC = ArrayContainer.none) {
	void update(AircraftStateT!AC ac_state, WakeT!AC wake, double dt) { assert(0); };
	Chunk inflow_at(immutable Vector!(4, Chunk) xyz) { assert(0); };
	void update_wing_circulation(WingStateT!AC wing_state) { assert(0); };
	void update_wing_dC_L(WingStateT!AC wing_state) { assert(0); };
	InducedVelocities compute_wing_induced_vel_on_blade(immutable Chunk x, immutable Chunk y, immutable Chunk z) { assert(0); };
	double wake_skew() { assert(0); };
	@nogc Frame* frame() { assert(0); };
}

class NullInflow(ArrayContainer AC) : InflowT!AC {

	alias RG = RotorGeometryT!AC;
	alias RS = RotorStateT!AC;
	alias RIS = RotorInputStateT!AC;

	package RG* rotor;
	package RIS* rotor_input;

	Frame* local_frame;

	double advance_ratio;
	double axial_advance_ratio;

	this(RG* _rotor, RIS* _rotor_input) {
		rotor = _rotor;
		rotor_input = _rotor_input;
		
		_rotor.frame.parent.children ~= new Frame(Vec3(1, 0, 0), PI, Vec3(0, 0, 0.0), _rotor.frame.parent, _rotor.frame.parent.name ~ " inflow", "connection");
		local_frame = _rotor.frame.parent.children[$-1];
	}

	override void update(AircraftStateT!AC ac_state, WakeT!AC wake, double dt) {
		auto inflow_local_freestream = local_frame.inverse_global_matrix * ac_state.freestream;
		
		advance_ratio = abs(inflow_local_freestream[0])/abs(rotor_input.angular_velocity*rotor.radius);
		axial_advance_ratio = inflow_local_freestream[2]/abs(rotor_input.angular_velocity*rotor.radius);

		auto global_inverse = local_frame.inverse_global_matrix;
		if (advance_ratio > 0) {
			immutable local_freestream = global_inverse*ac_state.freestream;
			immutable normal = Vec4(0, 0, -1, 0);
			immutable projected_freestream = local_freestream - local_freestream.dot(normal)*normal;
			immutable x_axis = Vec4(-1, 0, 0, 0);
			immutable double freestream_rotation = acos(projected_freestream.dot(x_axis)/projected_freestream.magnitude);
			local_frame.rotate(Vec3(0, 0, -1), freestream_rotation);
			local_frame.update(local_frame.parent.global_matrix);
		}

		inflow_local_freestream = local_frame.inverse_global_matrix * ac_state.freestream;
		
		advance_ratio = abs(inflow_local_freestream[0])/abs(rotor_input.angular_velocity*rotor.radius);
		axial_advance_ratio = inflow_local_freestream[2]/abs(rotor_input.angular_velocity*rotor.radius);
	}

	override Chunk inflow_at(immutable Vector!(4, Chunk) xyz) {
		immutable Chunk lambda = 0;
		return lambda;
	}

	override void update_wing_circulation(WingStateT!AC wing_state) {}
	override void update_wing_dC_L(WingStateT!AC wing_state) {}
	
	override double wake_skew() {
		return atan2(advance_ratio, axial_advance_ratio);
	}

	override InducedVelocities compute_wing_induced_vel_on_blade(immutable Chunk x, immutable Chunk y, immutable Chunk z) {
		InducedVelocities vel;
		vel.v_x[] = 0;
		vel.v_y[] = 0;
		vel.v_z[] = 0;
		return vel;
	}
	
	override @nogc Frame* frame() { return local_frame; }
}
