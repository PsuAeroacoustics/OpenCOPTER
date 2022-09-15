module opencopter.inflow.beddos;

import opencopter.aircraft;
import opencopter.inflow;
import opencopter.math;
import opencopter.memory;
import opencopter.wake;

import std.array;
import std.algorithm;
import std.math;
import std.range;

private auto inflow_at_impl(C)(BeddosInflow bi, auto ref C x, auto ref C y) {
	version(LDC) pragma(inline, true);
	static import std.math;
	import std.math : PI;

	Chunk inflow;
	
	immutable Chunk y_abs = abs(y);

	immutable Chunk x_minus_y_cube = -x[] - y_abs[]*y_abs[]*y_abs[];

	inflow[] = bi.k0 + bi.kx*x_minus_y_cube[] + bi.ky*y[];
	inflow[] *= bi.ffi;

	return inflow;
}

class BeddosInflow : Inflow {

	double ffi;
	private double kx;
	private double ky;
	private double k0;
	private double xi;
	//private double rotation_dir = 1.0;

	this(/+Direction direction+/) {
		/+if(direction == Direction.clockwise) {
			rotation_dir = 1.0;
		} else {
			rotation_dir = -1.0;
		}+/
	}

	@nogc double get_average_inflow() {
		return ffi;
	}

	@nogc void update(double C_T, ref RotorInputState rotor, ref RotorState rotor_state, double advance_ratio, double axial_advance_ratio, double dt) {
		update_impl(C_T, rotor, advance_ratio, axial_advance_ratio, dt);
	}

	@nogc void update(double C_T, RotorInputState* rotor, RotorState* rotor_state, double advance_ratio, double axial_advance_ratio, double dt) {
		update_impl(C_T, rotor, advance_ratio, axial_advance_ratio, dt);
	}

	private @nogc void update_impl(RIS)(double C_T, auto ref RIS rotor, double advance_ratio, double axial_advance_ratio, double dt)
		if(is_rotor_input_state!RIS)
	{
		import std.stdio : writeln;

		static import std.math;
		import std.math : sqrt, PI;
		import std.numeric : findRoot;

		double hover_inflow = sqrt(0.5*C_T);

		double forward_flight_inflow(double _ffi) {
			
			auto axial_plus_ffi = axial_advance_ratio + _ffi;
			auto ret = _ffi - hover_inflow*hover_inflow/(sqrt(advance_ratio*advance_ratio + axial_plus_ffi*axial_plus_ffi));
			
			return ret;
		}

		ffi = findRoot(&forward_flight_inflow, -0.7, 0.7);

		xi = atan2(advance_ratio, (axial_advance_ratio + ffi));

		kx = xi/2;
		k0 = 1.0 + 8.0*kx/(PI*15.0);
		ky = 0;
	}

	@nogc double wake_skew() {
		return xi;
	}
	
	@nogc Chunk inflow_at(immutable Chunk x, immutable Chunk y, immutable Chunk z, immutable Chunk x_e, double angle_of_attack) {

		Chunk inflow = this.inflow_at_impl(x, y)[];
		
		foreach(c_idx; 0..Chunk.length) {
			if(abs(y[c_idx]) > 1.01 || x[c_idx] > 1.0) {
				inflow[c_idx] = 0;
			} else if((abs(x[c_idx]) > x_e[c_idx])) {
				immutable y_abs = abs(y[c_idx]);

				immutable y_cube = y_abs*y_abs*y_abs;

				inflow[c_idx] = k0 - kx*y_cube + ky*y[c_idx];
				inflow[c_idx] *= ffi;
			}
		}

		return inflow;		
	}

	@nogc Chunk inflow_at(immutable Chunk r, immutable double cos_azimuth, immutable double sin_azimuth) {
		static import std.math;
		import std.math : PI;
	
		immutable Chunk x = r[]*cos_azimuth;
		immutable Chunk y = r[]*sin_azimuth;
		return this.inflow_at_impl(x, y);
	}
}
