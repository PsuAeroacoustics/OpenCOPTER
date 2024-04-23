module opencopter.inflow;

import opencopter.aircraft;
import opencopter.memory;

public import opencopter.inflow.beddos;
public import opencopter.inflow.huangpeters;

enum Direction {
	clockwise,
	counter_clockwise
}

interface Inflow {
	@nogc void update(double C_T, ref RotorInputState rotor, ref RotorState rotor_state, double advance_ratio, double axial_advance_ratio, double dt);
	@nogc void update(double C_T, RotorInputState* rotor, RotorState* rotor_state, double advance_ratio, double axial_advance_ratio, double dt);
	@nogc Chunk inflow_at(immutable Chunk x, immutable Chunk y, immutable Chunk z, immutable Chunk x_e, double angle_of_attack);
	@nogc Chunk inflow_at(immutable Chunk r, immutable double cos_azimuth, immutable double sin_azimuth);
	@nogc double wake_skew();
	//@nogc Frame* frame();
}
