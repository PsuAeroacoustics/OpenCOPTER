module opencopter.inflow.huangpeters.system;

import opencopter.aircraft;
import opencopter.config;
import opencopter.inflow;
import opencopter.inflow.huangpeters;
import opencopter.inflow.huangpeters.iteration;
import opencopter.inflow.huangpeters.math;
import opencopter.inflow.huangpeters.velocities;
import opencopter.math;
import opencopter.math.blas;
import opencopter.math.lapacke;
import opencopter.memory;
import opencopter.wake;

import numd.calculus.integration.forwardeuler;
import numd.calculus.integration.rk4;
import numd.linearalgebra.matrix;
import numd.utility : linspace;

import std.array;
import std.algorithm;
import std.complex;
import std.conv;
import std.math;
import std.range;
import std.stdio;
import std.typecons;

void build_vlm_matrix(ArrayContainer AC)(HuangPetersInflowT!AC infl, double advance_ratio, double axial_advance_ratio) {
	
	if(infl.average_inflow != 0) {
		infl.chi = atan2(advance_ratio, (infl.average_inflow + axial_advance_ratio));
		//debug writeln("advance_ratio: ", advance_ratio, " infl.average_inflow: ", infl.average_inflow, " axial_advance_ratio: ", axial_advance_ratio, " infl.chi: ", infl.chi*(180.0/PI));
	}

	infl.sin_chi = sin(infl.chi);
	infl.cos_chi = cos(infl.chi);
	infl.tan_chi = tan(infl.chi);

	immutable X = tan(-infl.chi/2.0);

	iterate_whole_matrix!(
		// odd-odd
		(r, j, m, n, row_idx, col_idx) {
			immutable l = min(r, m);
			immutable exp1 = abs(r - m);
			immutable exp2 = abs(r + m);

			if(r == 0) {
				infl.L_c[row_idx][col_idx] = X^^(m.to!double)*infl.Gamma[row_idx][col_idx];
			} else {
				infl.L_c[row_idx][col_idx] = (X^^exp1 + ((-1.0)^^l)*X^^exp2)*infl.Gamma[row_idx][col_idx];
			}
		},
		// odd-even
		(r, j, m, n, row_idx, col_idx) {
			immutable l = min(r, m);
			immutable exp1 = abs(r - m);
			immutable exp2 = abs(r + m);

			if(r == 0) {
				infl.L_c[row_idx][col_idx] = X^^(m.to!double)*infl.Gamma[row_idx][col_idx];
			} else {
				infl.L_c[row_idx][col_idx] = (X^^exp1 + ((-1.0)^^l)*X^^exp2)*infl.Gamma[row_idx][col_idx];
			}
		},
		// even-odd
		(r, j, m, n, row_idx, col_idx) {
			immutable l = min(r, m);
			immutable exp1 = abs(r - m);
			immutable exp2 = abs(r + m);

			if(r == 0) {
				infl.L_c[row_idx][col_idx] = X^^m*infl.Gamma[row_idx][col_idx];
			} else {
				infl.L_c[row_idx][col_idx] = (X^^exp1 + ((-1.0)^^l)*X^^exp2)*infl.Gamma[row_idx][col_idx];
			}
		},
		// even-even
		(r, j, m, n, row_idx, col_idx) {
			immutable l = min(r, m);
			immutable exp1 = abs(r - m);
			immutable exp2 = abs(r + m);

			if(r == 0) {
				infl.L_c[row_idx][col_idx] = X^^m*infl.Gamma[row_idx][col_idx];
			} else {
				infl.L_c[row_idx][col_idx] = (X^^exp1 + ((-1.0)^^l)*X^^exp2)*infl.Gamma[row_idx][col_idx];
			}
		}
	)(infl.Mo, infl.Me);


	iterate_whole_matrix_sin!(
		// odd-odd
		(r, j, m, n, row_idx, col_idx) {

			immutable double l = min(r, m);
			immutable double exp1 = abs(r - m);
			immutable double exp2 = abs(r + m);

			immutable gamma_col_idx = col_idx + odd_states[infl.Mo][0].length;
			immutable gamma_row_idx = row_idx + odd_states[infl.Mo][0].length;

			infl.L_s[row_idx][col_idx] = (X^^exp1 - ((-1.0)^^l)*X^^exp2)*infl.Gamma[gamma_row_idx][gamma_col_idx];
		},
		// odd-even
		(r, j, m, n, row_idx, col_idx) {
			if(r == 0) {
				infl.L_s[row_idx][col_idx] = 0;
			} else {
				immutable double l = min(r, m);
				immutable double exp1 = abs(r - m);
				immutable double exp2 = abs(r + m);

				immutable gamma_col_idx = col_idx + infl.total_odd_states + even_states[infl.Me][0].length - infl.total_odd_sin_states;
				immutable gamma_row_idx = row_idx + odd_states[infl.Mo][0].length;

				infl.L_s[row_idx][col_idx] = (X^^exp1 - ((-1.0)^^l)*X^^exp2)*infl.Gamma[gamma_row_idx][gamma_col_idx];
			}
		},
		// even-odd
		(r, j, m, n, row_idx, col_idx) {

			if(r == 0) {

			} else {
				immutable l = min(r, m);
				immutable exp1 = abs(r - m);
				immutable exp2 = abs(r + m);

				immutable gamma_col_idx = col_idx + odd_states[infl.Mo][0].length;
				immutable gamma_row_idx = row_idx + infl.total_odd_states + even_states[infl.Me][0].length - infl.total_odd_sin_states;

				infl.L_s[row_idx][col_idx] = (X^^exp1 - ((-1.0)^^l)*X^^exp2)*infl.Gamma[gamma_row_idx][gamma_col_idx];
			}
		},
		// even-even
		(r, j, m, n, row_idx, col_idx) {
			if(r == 0) {
				infl.L_s[row_idx][col_idx] = 0;
			} else {
				immutable l = min(r, m);
				immutable exp1 = abs(r - m);
				immutable exp2 = abs(r + m);

				immutable gamma_col_idx = col_idx + infl.total_odd_states + even_states[infl.Me][0].length - infl.total_odd_sin_states;
				immutable gamma_row_idx = row_idx + infl.total_odd_states + even_states[infl.Me][0].length - infl.total_odd_sin_states;

				infl.L_s[row_idx][col_idx] = (X^^exp1 - ((-1.0)^^l)*X^^exp2)*infl.Gamma[gamma_row_idx][gamma_col_idx];
			}
		}
	)(infl.Mo, infl.Me);

	foreach(r_idx, ref L_c_row; infl.L_c) {
		infl.L_c_inv[r_idx][] = L_c_row[];
	}

	foreach(r_idx, ref L_s_row; infl.L_s) {
		infl.L_s_inv[r_idx][] = L_s_row[];
	}

	// Invert L_c matrix
	int info = 0;
	info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, infl.total_states, infl.total_states, infl.L_c_inv[0].ptr, infl.total_states, infl.ipiv.ptr);
	assert(info == 0, "Failed to invert L_c matrix");
	info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, infl.total_states, infl.L_c_inv[0].ptr, infl.total_states, infl.ipiv.ptr);
	assert(info == 0, "Failed to invert L_c matrix");

	info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, infl.total_sin_states, infl.total_sin_states, infl.L_s_inv[0].ptr, infl.total_sin_states, infl.ipiv.ptr);
	assert(info == 0, "Failed to invert L_s matrix");
	info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, infl.total_sin_states, infl.L_s_inv[0].ptr, infl.total_sin_states, infl.ipiv.ptr);
	assert(info == 0, "Failed to invert L_s matrix");

	// Multiply L_c_inv*M_c
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, infl.total_states, infl.total_states, infl.total_states, 1.0, infl.L_c_inv[0].ptr, infl.total_states, infl.M_c[0].ptr, infl.total_states, 0.0, infl.VLM_c[0].ptr, infl.total_states);

	// Multiply L_s_inv*M_s
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, infl.total_sin_states, infl.total_sin_states, infl.total_sin_states, 1.0, infl.L_s_inv[0].ptr, infl.total_sin_states, infl.M_s[0].ptr, infl.total_sin_states, 0.0, infl.VLM_s[0].ptr, infl.total_sin_states);

	

	immutable v_inf_sin_squared = (advance_ratio*infl.sin_chi)^^2.0;
	immutable v_inf_cos = advance_ratio*infl.cos_chi;
	immutable double Vt = sqrt(v_inf_sin_squared + (v_inf_cos + infl.average_inflow)^^2.0);
	immutable double V = (v_inf_sin_squared + (v_inf_cos + infl.average_inflow)*(v_inf_cos + 2.0*infl.average_inflow))/Vt;

	infl.LM_c_scratch[] = infl.VLM_c[0][];// * alpha[];

	// Build the non linearity matrix
	infl.VLM_c[0][] *= Vt;

	foreach(idx, ref vlm; infl.VLM_c[1..$]) {
		vlm[] *= V;
	}

	foreach(idx, ref vlm; infl.VLM_s[0..$]) {
		vlm[] *= V;
	}

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, infl.total_states, infl.total_states, infl.total_states, 1.0, infl.M_c_inv[0].ptr, infl.total_states, infl.L_c[0].ptr, infl.total_states, 0.0, infl.VLM_c_inv[0].ptr, infl.total_states);
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, infl.total_sin_states, infl.total_sin_states, infl.total_sin_states, 1.0, infl.M_s_inv[0].ptr, infl.total_sin_states, infl.L_s[0].ptr, infl.total_sin_states, 0.0, infl.VLM_s_inv[0].ptr, infl.total_sin_states);

	foreach(idx, ref vlM_c_inv; infl.VLM_c_inv) {
		vlM_c_inv[0] *= 1.0/Vt;
		vlM_c_inv[1..$] *= 1.0/V;
	}

	foreach(idx, ref vlM_s_inv; infl.VLM_s_inv) {
		vlM_s_inv[0..$] *= 1.0/V;
	}

	foreach(idx, ref vlM_c_inv; infl.VLM_c_inv) {
		vlM_c_inv[] *= infl.adjoint_mat[];
	}

	foreach(idx, ref vlM_s_inv; infl.VLM_s_inv) {
		vlM_s_inv[] *= infl.adjoint_mat_sin[];
	}

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, infl.total_states, infl.total_states, infl.total_states, 1.0, infl.VLM_c_inv[0].ptr, infl.total_states, infl.VLM_c[0].ptr, infl.total_states, 0.0, infl.QS_c_mat[0].ptr, infl.total_states);

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, infl.total_sin_states, infl.total_sin_states, infl.total_sin_states, 1.0, infl.VLM_s_inv[0].ptr, infl.total_sin_states, infl.VLM_s[0].ptr, infl.total_sin_states, 0.0, infl.QS_s_mat[0].ptr, infl.total_sin_states);

}

void simple_harmonic_solution(ArrayContainer AC)(HuangPetersInflowT!AC infl, double advance_ratio, double axial_advance_ratio) {

	infl.build_vlm_matrix(advance_ratio, axial_advance_ratio);

	int info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, infl.total_states, infl.total_states, infl.VLM_c[0].ptr, infl.total_states, infl.ipiv.ptr);
	assert(info == 0, "Failed to invert VLM_c matrix");
	info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, infl.total_states, infl.VLM_c[0].ptr, infl.total_states, infl.ipiv.ptr);
	assert(info == 0, "Failed to invert VLM_c matrix");


	info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, infl.total_sin_states, infl.total_sin_states, infl.VLM_s[0].ptr, infl.total_sin_states, infl.ipiv.ptr);
	assert(info == 0, "Failed to invert VLM_s matrix");
	info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, infl.total_sin_states, infl.VLM_s[0].ptr, infl.total_sin_states, infl.ipiv.ptr);
	assert(info == 0, "Failed to invert VLM_s matrix");


	cblas_dgemv(CblasRowMajor, CblasNoTrans, infl.total_states, infl.total_states, 1.0, infl.VLM_c[0].ptr, infl.total_states, infl.tau_c.ptr, 1, 0.0, infl.state_history[infl.get_circular_index(infl.curr_state)][0..infl.total_states].ptr, 1);

	cblas_dgemv(CblasRowMajor, CblasNoTrans, infl.total_sin_states, infl.total_sin_states, 1.0, infl.VLM_c[0].ptr, infl.total_sin_states, infl.tau_s.ptr, 1, 0.0, infl.state_history[infl.get_circular_index(infl.curr_state)][2*infl.total_states..2*infl.total_states + infl.total_sin_states].ptr, 1);
	
	infl.alpha = infl.state_history[infl.get_circular_index(infl.curr_state)][0..infl.total_states];
	infl.beta = infl.state_history[infl.get_circular_index(infl.curr_state)][2*infl.total_states..2*infl.total_states + infl.total_sin_states];

	foreach(i; 0..infl.total_states) {
		infl.state_history[infl.get_circular_index(infl.curr_state)][infl.total_states + i] = 0;
		foreach(k; 0..infl.total_states) {
			infl.state_history[infl.get_circular_index(infl.curr_state)][infl.total_states + i] += infl.QS_c_mat[i][k]*infl.alpha[k];
		}
	}

	foreach(i; 0..infl.total_sin_states) {
		infl.state_history[infl.get_circular_index(infl.curr_state)][2*infl.total_states + infl.total_sin_states + i] = 0;
		foreach(k; 0..infl.total_sin_states) {
			infl.state_history[infl.get_circular_index(infl.curr_state)][2*infl.total_states + infl.total_sin_states + i] += infl.QS_s_mat[i][k]*infl.beta[k];
		}
	}

	infl.lambda = infl.state_history[infl.get_circular_index(infl.curr_state)][infl.total_states..2*infl.total_states];

	infl.lambda_s = infl.state_history[infl.get_circular_index(infl.curr_state)][2*infl.total_states + infl.total_sin_states..$];

	infl.a[] = 0;
	infl.b[] = 0;
	// Get MD state variables
	cblas_dgemv(CblasRowMajor, CblasNoTrans, infl.total_odd_states, infl.total_odd_states, 1.0, infl.A_inv[0].ptr, infl.total_odd_states, infl.alpha.ptr, 1, 0.0, infl.a.ptr, 1);

	cblas_dgemv(CblasRowMajor, CblasNoTrans, infl.total_odd_sin_states, infl.total_odd_sin_states, 1.0, infl.A_s_inv[0].ptr, infl.total_odd_sin_states, infl.beta.ptr, 1, 0.0, infl.b.ptr, 1);

	// Get MD adjoint variables
	cblas_dgemv(CblasRowMajor, CblasNoTrans, infl.total_odd_states, infl.total_odd_states, 1.0, infl.A_inv[0].ptr, infl.total_odd_states, infl.lambda.ptr, 1, 0.0, infl.delta.ptr, 1);

	cblas_dgemv(CblasRowMajor, CblasNoTrans, infl.total_odd_sin_states, infl.total_odd_sin_states, 1.0, infl.A_s_inv[0].ptr, infl.total_odd_sin_states, infl.lambda_s.ptr, 1, 0.0, infl.delta_s.ptr, 1);

	infl.a[infl.total_odd_states..$] = infl.alpha[infl.total_odd_states..$];
	infl.delta[infl.total_odd_states..$] = infl.lambda[infl.total_odd_states..$];

	infl.b[infl.total_odd_sin_states..$] = infl.beta[infl.total_odd_sin_states..$];
	infl.delta_s[infl.total_odd_sin_states..$] = infl.lambda_s[infl.total_odd_sin_states..$];

}

void system_derivative(ArrayContainer AC, RIS, RS)(double[] state_dot, double[] state, double t, double dt, HuangPetersInflowT!AC infl, auto ref RIS rotor, auto ref RS rotor_state, double advance_ratio, double axial_advance_ratio) {
	
	double[] alpha = state[0..infl.total_states];
	double[] lambda = state[infl.total_states..2*infl.total_states];

	double[] beta = state[2*infl.total_states..2*infl.total_states + infl.total_sin_states];
	double[] lambda_s = state[2*infl.total_states + infl.total_sin_states..$];

	build_vlm_matrix(infl, advance_ratio, axial_advance_ratio);

	// scratch * alpha = V*L_c_inv*M_c*alpha
	cblas_dgemv(CblasRowMajor, CblasNoTrans, infl.total_states, infl.total_states, 1.0, infl.VLM_c[0].ptr, infl.total_states, alpha.ptr, 1, 0.0, infl.alpha_scratch.ptr, 1);

	// scratch * beta = V*L_s_inv*M_s*beta
	cblas_dgemv(CblasRowMajor, CblasNoTrans, infl.total_sin_states, infl.total_sin_states, 1.0, infl.VLM_s[0].ptr, infl.total_sin_states, beta.ptr, 1, 0.0, infl.beta_scratch.ptr, 1);

	// scratch * delta = V*L_c_inv*M_c*delta
	cblas_dgemv(CblasRowMajor, CblasNoTrans, infl.total_states, infl.total_states, 1.0, infl.VLM_c[0].ptr, infl.total_states, lambda.ptr, 1, 0.0, infl.delta_scratch.ptr, 1);

	// scratch * delta_s = V*L_s_inv*M_s*delta_s
	cblas_dgemv(CblasRowMajor, CblasNoTrans, infl.total_sin_states, infl.total_sin_states, 1.0, infl.VLM_s[0].ptr, infl.total_sin_states, lambda_s.ptr, 1, 0.0, infl.delta_s_scratch.ptr, 1);

	infl.alpha_scratch[] = infl.tau_c[] - infl.alpha_scratch[];
	infl.beta_scratch[] = infl.tau_s[] - infl.beta_scratch[];

	// M_c_inv*D*alpha_scratch = M_c_inv*D*(tau_c - V*L_c_inv*M_c*alpha)
	cblas_dgemv(CblasRowMajor, CblasNoTrans, infl.total_states, infl.total_states, 1.0, infl.M_c_inv_D[0].ptr, infl.total_states, infl.alpha_scratch.ptr, 1, 0.0, state_dot[0..infl.total_states].ptr, 1);

	cblas_dgemv(CblasRowMajor, CblasNoTrans, infl.total_sin_states, infl.total_sin_states, 1.0, infl.M_s_inv_D[0].ptr, infl.total_sin_states, infl.beta_scratch.ptr, 1, 0.0, state_dot[2*infl.total_states..2*infl.total_states + infl.total_sin_states].ptr, 1);
}

class HuangPetersInflowT(ArrayContainer AC = ArrayContainer.none) : Inflow {

	alias RG = RotorGeometryT!AC;
	alias RS = RotorStateT!AC;
	alias RIS = RotorInputStateT!AC;

	package immutable long Mo;
	package immutable long Me;
	package long N;

	package Chunk[] mpsi_buff;
	package Chunk[] sin_mpsi_buff;

	package double[][] M_c; // Constant.
	package double[][] L_c_inv; // Modification on Gamma by skew angle then inverted each timestep.
	package double[][] L_c; // Modification on Gamma by skew angle then inverted each timestep.
	package double[][] M_s; // Constant.
	package double[][] L_s_inv; // Modification on Gamma by skew angle then inverted each timestep.
	package double[][] L_s; // Modification on Gamma by skew angle then inverted each timestep.
	package double[][] Gamma; // Constant throughout.
	package double[][] M_c_inv_D; // Damping. Constant
	package double[][] M_s_inv_D; // Damping. Constant
	package double[][] VLM_c; // This holds the result of the matrix mults each timestep.
	package double[][] VLM_c_inv; // This holds the result of the matrix mults each timestep.
	package double[][] VLM_s; // This holds the result of the matrix mults each timestep.
	package double[][] VLM_s_inv; // This holds the result of the matrix mults each timestep.
	double[][] QS_c_mat; // This holds the result of the matrix mults each timestep.
	double[][] QS_s_mat; // This holds the result of the matrix mults each timestep.
	package double[][] A_inv;
	package double[][] A;
	package double[][] A_s;
	package double[][] A_s_inv;
	package double[][] M_c_inv;
	package double[][] M_s_inv;
	package double[][] D;
	package double[][] D_s;

	package int[] ipiv; // Used by lapack for matrix inversion. Keep it around instead of re-allocing each update

	package ptrdiff_t z_length;
	package ptrdiff_t time_history;
	package ptrdiff_t curr_state;
	package Chunk[] time_delay_alpha;
	package Chunk[] time_delay_alpha_2;
	package Chunk[] time_delay_a;
	package Chunk[] time_delay_a_2;
	package Chunk[] time_delay_beta;
	package Chunk[] time_delay_beta_2;
	package Chunk[] time_delay_b;
	package Chunk[] time_delay_b_2;
	package double[][] state_history;
	package double[] times;
	double[] tau_c; // Pressure coefficients
	double[] tau_s; // Pressure coefficients
	package double[] tau_c_scratch;
	package double[] tau_s_scratch;
	double[] alpha; // NH state variables
	double[] beta; // NH state variables
	package double[] a; // MD state variables
	package double[] b; // MD state variables
	package double[] alpha_scratch;
	package double[] beta_scratch;
	double[] delta; // MD adjoint variables
	double[] lambda; // NH adjoint variables
	double[] delta_s; // MD adjoint variables
	double[] lambda_s; // NH adjoint variables
	package double[] delta_scratch;
	package double[] delta_s_scratch;
	package double[] adjoint_mat;
	package double[] adjoint_mat_sin;
	package double[] LM_c_scratch;
	package double[][] K_table;
	package double[] average_inflow_array;

	int total_states;
	int total_odd_states;
	int total_even_states;

	int total_sin_states;
	int total_odd_sin_states;
	int total_even_sin_states;

	package double average_inflow;

	package double[][] P_coefficients;
	package double[][] P_coefficients_nh;
	package Chunk[][] Qmn_bar;
	package Chunk[][] P_md_mn_bar;
	package Chunk[][] P_nh_mn_bar;
	
	alias Integrator = ForwardEuler!double;
	Integrator integrator;

	package Chunk[] blade_scratch;
	package Chunk[] blade_scratch_s;

	package RG* rotor;
	package RS* rotor_state;
	package RIS* rotor_input;

	package double advance_ratio;
	package double axial_advance_ratio;
	package double omega;

	package double chi;
	package double sin_chi;
	package double cos_chi;
	package double tan_chi;
	package double cot_chi;

	size_t n_r = 16;
	size_t n_psi = 6;
	double v_0;

	Chunk[] contraction_array;
	double[] contraction_array_alias;
	Chunk[] contraction_z_array;
	double[] contraction_z_array_alias;

	bool contraction_mapping = false;
	bool debug_coords = false;

	double[] zero;
	Chunk[] c_zero;

	Frame* local_frame;
	Mat4 global_inverse;

	@nogc Frame* frame() {
		return local_frame;
	}

	@nogc Mat4 inverse_global_frame() {
		return global_inverse;
	}

	@nogc package ptrdiff_t get_circular_index(ptrdiff_t idx) {
		return ((idx % time_history) + time_history) % time_history;
	}

	this(long _Mo, long _Me, RG* _rotor, RS* _rotor_state,RIS* _rotor_input, double dt) {

		_rotor.frame.parent.children ~= new Frame(Vec3(1, 0, 0), PI, /+_rotor.frame.local_position()+/Vec3(0, 0, 0.0), _rotor.frame.parent, _rotor.frame.parent.name ~ " inflow", "connection");
		local_frame = _rotor.frame.parent.children[$-1];
		local_frame.local_matrix[1, 1] *= -1.0;

		size_t len = round(2.0*PI/(dt*235.325)).to!size_t;
		//ai_idx = 0;
		average_inflow_array = new double[len];
		average_inflow_array[] = 0;
		rotor = _rotor;
		rotor_state = _rotor_state;
		rotor_input = _rotor_input;
		size_t num_chunks = rotor.blades[0].chunks.length;

		openblas_set_num_threads(1);

		z_length = 48;
		contraction_z_array_alias = linspace(-5.0, 5.0, z_length);
		contraction_z_array = cast(Chunk[])contraction_z_array_alias.ptr[0..z_length];

		contraction_array_alias = new double[z_length];
		contraction_array = cast(Chunk[])contraction_array_alias.ptr[0..z_length];

		blade_scratch = new Chunk[num_chunks];
		blade_scratch_s = new Chunk[num_chunks];

		Mo = _Mo;
		Me = _Me;

		mpsi_buff = new Chunk[max(Mo, Me) + 1];
		sin_mpsi_buff = new Chunk[max(Mo, Me) + 1];

		N = 0;

		total_states = 0;
		total_even_states = 0;
		total_odd_states = 0;

		iterate_odds!((m, n, idx) {

			size_t num_coefficients = n - m + 2;
			P_coefficients ~= new double[num_coefficients];
			P_coefficients[$-1][0] = (-1.0)^^m.to!double*((-1.0)^^m.to!double)*(2.0^^n.to!double);
			foreach(k; m..n + 1) {
				P_coefficients[$-1][k - m + 1] = (k.factorial/(k - m).factorial)*(n.falling_factorial(k)/k.factorial)*((0.5*(n.to!double + k.to!double - 1.0)).falling_factorial(n)/n.factorial );
			}

			N = max(n + 1, N);
			immutable rho = sqrt(1.0/(2.0*n.to!double + 1.0)*(n + m).factorial/(n - m).factorial);
			P_coefficients[$-1][0] /= rho;

			//writeln("idx: ", idx, "; m: ", m, ", n: ", n);
			total_states++;
			total_odd_states++;
		})(Mo, 0);

		iterate_evens!((m, n, idx) {
			size_t num_coefficients = n - m + 2;
			P_coefficients ~= new double[num_coefficients];
			P_coefficients[$-1][0] = (-1.0)^^m.to!double*((-1.0)^^m.to!double)*(2.0^^n.to!double);
			foreach(k; m..n + 1) {
				P_coefficients[$-1][k - m + 1] = (k.factorial/(k - m).factorial)*(n.falling_factorial(k)/k.factorial)*((0.5*(n.to!double + k.to!double - 1.0)).falling_factorial(n)/n.factorial );
			}

			N = max(n + 1, N);

			immutable rho = sqrt(1.0/(2.0*n.to!double + 1.0)*(n + m).factorial/(n - m).factorial);
			P_coefficients[$-1][0] /= rho;

			total_states++;
			total_even_states++;
		})(Me, total_odd_states);

		iterate_odds!((r, j, idx) {

			size_t num_coefficients = j - r + 2;
			P_coefficients_nh ~= new double[num_coefficients];
			P_coefficients_nh[$-1][0] = sqrt((2.0*j.to!double + 1.0)*H(r, j));

			foreach(qidx, q; iota(r, j, 2).enumerate) {
				P_coefficients_nh[$-1][q - r + 1] = (-1.0)^^(0.5*(q.to!double - r.to!double))*double_factorial(j + q);
				P_coefficients_nh[$-1][q - r + 1] /= double_factorial(q - r)*double_factorial(q + r)*double_factorial(j - q - 1);
			}
		})(Mo, 0);

		iterate_evens!((r, j, idx) {

			size_t num_coefficients = j - r + 2;
			P_coefficients_nh ~= new double[num_coefficients];
			P_coefficients_nh[$-1][0] = sqrt((2.0*j.to!double + 1.0)*H(r, j));

			foreach(qidx, q; iota(r, j, 2).enumerate) {
				P_coefficients_nh[$-1][q - r + 1] = (-1.0)^^(0.5*(q.to!double - r.to!double))*double_factorial(j + q);
				P_coefficients_nh[$-1][q - r + 1] /= double_factorial(q - r)*double_factorial(q + r)*double_factorial(j - q - 1);
			}
		})(Me, 0);

		total_even_sin_states = total_even_states - even_states[Me][0].length.to!int;
		total_odd_sin_states = total_odd_states - odd_states[Mo][0].length.to!int;

		total_sin_states = total_odd_sin_states + total_even_sin_states;

		debug writeln("total_states: ", total_states);
		debug writeln("total_odd_states: ", total_odd_states);
		debug writeln("total_even_states: ", total_even_states);
		debug writeln("total_sin_states: ", total_sin_states);
		debug writeln("total_odd_sin_states: ", total_odd_sin_states);
		debug writeln("total_even_sin_states: ", total_even_sin_states);
		
		integrator = new Integrator(2*total_states + 2*total_sin_states);

		double[][] Meo = allocate_dense(total_even_states, total_odd_states);
		double[][] Meo_new = allocate_dense(total_even_states, total_odd_states);
		double[][] Meo_s = allocate_dense(total_even_sin_states, total_odd_sin_states);
		double[][] Meo_s_new = allocate_dense(total_even_sin_states, total_odd_sin_states);
		K_table = allocate_dense(max(Me, Mo) + 1, N);
		D = allocate_dense(total_states, total_states);
		D_s = allocate_dense(total_sin_states, total_sin_states);
		QS_c_mat = allocate_dense(total_states, total_states);
		M_c_inv = allocate_dense(total_states, total_states);
		QS_s_mat = allocate_dense(total_sin_states, total_sin_states);
		M_s_inv = allocate_dense(total_sin_states, total_sin_states);
		A_inv = allocate_dense(total_odd_states, total_odd_states);
		A_s_inv = allocate_dense(total_odd_sin_states, total_odd_sin_states);
		A_s = allocate_dense(total_odd_sin_states, total_odd_sin_states);
		A = allocate_dense(total_odd_states, total_odd_states);
		VLM_c = allocate_dense(total_states, total_states);
		VLM_c_inv = allocate_dense(total_states, total_states);
		M_c = allocate_dense(total_states, total_states);
		M_c_inv_D = allocate_dense(total_states, total_states);
		L_c_inv = allocate_dense(total_states, total_states);
		L_c = allocate_dense(total_states, total_states);
		VLM_s = allocate_dense(total_sin_states, total_sin_states);
		VLM_s_inv = allocate_dense(total_sin_states, total_sin_states);
		M_s = allocate_dense(total_sin_states, total_sin_states);
		M_s_inv_D = allocate_dense(total_sin_states, total_sin_states);
		L_s_inv = allocate_dense(total_sin_states, total_sin_states);
		L_s = allocate_dense(total_sin_states, total_sin_states);
		Gamma = allocate_dense(total_states, total_states);

		foreach(ref K; K_table) {
			K[] = 0;
		}

		time_history = 1_000_000;
		state_history = allocate_dense(time_history, 2*total_states + 2*total_sin_states);
		
		Qmn_bar = allocate_dense_chunk(max(Me, Mo) + 1, N);
		P_nh_mn_bar = allocate_dense_chunk(max(Me, Mo) + 1, N);
		P_md_mn_bar = allocate_dense_chunk(max(Me, Mo) + 1, N);

		times = new double[time_history];
		
		zero = new double[total_sin_states];
		c_zero = new Chunk[total_sin_states];
		time_delay_alpha = new Chunk[2*total_states];
		time_delay_alpha_2 = new Chunk[2*total_states];
		time_delay_a = new Chunk[2*total_states];
		time_delay_a_2 = new Chunk[2*total_states];

		time_delay_beta = new Chunk[2*total_sin_states];
		time_delay_beta_2 = new Chunk[2*total_sin_states];
		time_delay_b = new Chunk[2*total_sin_states];
		time_delay_b_2 = new Chunk[2*total_sin_states];

		a = new double[total_states];
		b = new double[total_sin_states];
		delta = new double[total_states];
		delta_s = new double[total_sin_states];
		alpha_scratch = new double[total_states];
		beta_scratch = new double[total_sin_states];
		delta_scratch = new double[total_states];
		delta_s_scratch = new double[total_sin_states];
		tau_c = new double[total_states];
		tau_c_scratch = new double[total_states];
		tau_s = new double[total_sin_states];
		tau_s_scratch = new double[total_sin_states];
		adjoint_mat = new double[total_states];
		adjoint_mat_sin = new double[total_sin_states];
		LM_c_scratch = new double[total_states];

		zero[] = 0;
		foreach(ref ch; c_zero) {
			ch[] = 0;
		}
		alpha_scratch[] = 0;
		beta_scratch[] = 0;
		delta_scratch[] = 0;
		delta_s_scratch[] = 0;
		tau_c[] = 0;
		tau_c_scratch[] = 0;
		tau_s[] = 0;
		tau_s_scratch[] = 0;
		adjoint_mat[] = 0;
		adjoint_mat_sin[] = 0;

		state_history.zero_matrix;

		alpha = state_history[0][0..total_states];
		lambda = state_history[0][total_states..2*total_states];

		beta = state_history[0][2*total_states..2*total_states + total_sin_states];
		lambda_s = state_history[0][2*total_states + total_sin_states..$];

		A_inv.zero_matrix;
		A_s_inv.zero_matrix;
		Meo.zero_matrix;
		Meo_new.zero_matrix;
		Meo_s.zero_matrix;
		Meo_s_new.zero_matrix;
		M_c_inv.zero_matrix;
		M_c.zero_matrix;
		M_s_inv.zero_matrix;
		M_s.zero_matrix;
		D.zero_matrix;
		D_s.zero_matrix;
		Gamma.zero_matrix;
		L_c_inv.zero_matrix;
		VLM_c.zero_matrix;
		VLM_c_inv.zero_matrix;
		QS_c_mat.zero_matrix;
		L_s_inv.zero_matrix;
		VLM_s.zero_matrix;
		VLM_s_inv.zero_matrix;
		QS_s_mat.zero_matrix;

		iterate_even_odd!(
			(r, j, m, n, row_idx, col_idx) {
				if(r == m) {
					if((j == n - 1) || (j == n + 1)) {
						
						Meo[row_idx][col_idx] = 1;
						Meo[row_idx][col_idx] /= sqrt(H(m, n)*H(r, j))*sqrt((2.0*n.to!double + 1.0)*(2.0*j.to!double + 1.0));
					}
				}
			}
		)(Mo, Me, 0, 0);

		iterate_odd_odd!((r, j, m, n, row_idx, col_idx) {
			if(r == m) {
				A_inv[row_idx][col_idx] = 2.0*(-1.0)^^(0.5*(n + j - 2*r).to!double);
				A_inv[row_idx][col_idx] *= sqrt((2.0*n.to!double + 1.0)*(2.0*j.to!double + 1.0));
				A_inv[row_idx][col_idx] /= sqrt(H(r, n)*H(r, j))*(n + j).to!double*(n.to!double + j.to!double + 2.0)*((n.to!double - j.to!double)^^2.0 - 1.0);
			}
		})(Mo, 0, 0);

		foreach(r_idx; 0..total_odd_states) {
			A[r_idx][] = A_inv[r_idx][];
		}

		iterate_odd_odd_sin!((r, j, m, n, row_idx, col_idx) {
			if(r == m) {
				A_s_inv[row_idx][col_idx] = 2.0*(-1.0)^^(0.5*(n + j - 2*r).to!double);
				A_s_inv[row_idx][col_idx] *= sqrt((2.0*n.to!double + 1.0)*(2.0*j.to!double + 1.0));
				A_s_inv[row_idx][col_idx] /= sqrt(H(r, n)*H(r, j))*(n + j).to!double*(n.to!double + j.to!double + 2.0)*((n.to!double - j.to!double)^^2.0 - 1.0);
			}
		})(Mo, 0, 0);

		foreach(r_idx; 0..total_odd_sin_states) {
			A_s[r_idx][] = A_s_inv[r_idx][];
		}

		int info = 0;
		ipiv = new int[total_odd_states];

		info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, total_odd_states, total_odd_states, A_inv[0].ptr, total_odd_states, ipiv.ptr);
		assert(info == 0, "Failed to invert A matrix");
		info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, total_odd_states, A_inv[0].ptr, total_odd_states, ipiv.ptr);
		assert(info == 0, "Failed to invert A matrix");

		info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, total_odd_sin_states, total_odd_sin_states, A_s_inv[0].ptr, total_odd_sin_states, ipiv.ptr);
		assert(info == 0, "Failed to invert A_s matrix");
		info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, total_odd_sin_states, A_s_inv[0].ptr, total_odd_sin_states, ipiv.ptr);
		assert(info == 0, "Failed to invert A_s matrix");

		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, total_even_states, total_odd_states, total_odd_states, 1.0, Meo[0].ptr, total_odd_states, A_inv[0].ptr, total_odd_states, 0.0, Meo_new[0].ptr, total_odd_states);

		iterate_whole_matrix!(
			// odd-odd
			(r, j, m, n, row_idx, col_idx) {
				if(row_idx == col_idx) {
					M_c[row_idx][col_idx] = 1;
				}

				if(r == m) {
					D[row_idx][col_idx] = kronecker_delta(j, n)/K(m, n);
				}

				if(((r + m) % 2 != 0) && ((j == (n - 1)) || (j == (n + 1)))) {
					Gamma[row_idx][col_idx] = sgn(r - m);
					Gamma[row_idx][col_idx] /= sqrt(K(m, n)*K(r, j))*sqrt((2*n + 1).to!double*(2*j + 1).to!double);
				}

				if((r + m) % 2 == 0) {
					Gamma[row_idx][col_idx] = 2.0*(-1.0)^^(0.5*(n + j - 2.0*r));
					Gamma[row_idx][col_idx] *= sqrt((2.0*n + 1.0)*(2.0*j + 1.0));
					Gamma[row_idx][col_idx] /= sqrt(H(m, n)*H(r, j))*(n + j).to!double*(n + j + 2).to!double*((n - j).to!double^^2.0 - 1.0);
				}
			},
			// odd-even
			(r, j, m, n, row_idx, col_idx) {
				if(r == m) {
					if((j == n - 1) || (j == n + 1)) {
						M_c[row_idx][col_idx] = 1;
						M_c[row_idx][col_idx] /= sqrt(H(m, n)*H(r, j))*sqrt((2*n + 1).to!double*(2*j + 1).to!double);
					}

					D[row_idx][col_idx] = 2.0*sqrt((2*j + 1).to!double*(2*n + 1).to!double);
					D[row_idx][col_idx] *= (-1.0)^^(0.5*(j + 3*n - 1).to!double);
					D[row_idx][col_idx] /= PI*sqrt(H(m, n)*H(m, j))*(j + n + 1).to!double*(j - n).to!double;
				}

				if((r + m) % 2 != 0) {
					Gamma[row_idx][col_idx] = sgn(r - m)*4.0*(-1.0)^^(0.5*(3.0*n + j + 2.0*m - 2.0*r));
					Gamma[row_idx][col_idx] *= sqrt((2.0*n + 1.0)*(2.0*j + 1.0));
					Gamma[row_idx][col_idx] /= PI*sqrt(H(m, n)*H(r, j))*(n + j)*(n + j + 2.0)*((n - j)^^2.0 - 1.0);
				}

				if(((r + m) % 2 == 0) && ((j == (n - 1)) || (j == (n + 1)))) {
					Gamma[row_idx][col_idx] = 1.0;
					Gamma[row_idx][col_idx] /= sqrt(H(m, n)*H(r, j))*sqrt((2*n + 1).to!double*(2*j + 1).to!double);
				}
			},
			// even-odd
			(r, j, m, n, row_idx, col_idx) {
				if(r == m) {
					D[row_idx][col_idx] = 2.0*sqrt((2.0*j.to!double + 1)*(2.0*n.to!double + 1));
					D[row_idx][col_idx] *= (-1.0)^^(0.5*(j + 3*n - 1).to!double);
					D[row_idx][col_idx] /= PI*sqrt(H(m, n)*H(m, j))*(j + n + 1).to!double*(j - n).to!double;
				}

				M_c[row_idx][col_idx] = Meo_new[row_idx - total_odd_states][col_idx];

				if((r + m) % 2 != 0) {
					Gamma[row_idx][col_idx] = sgn(r - m)*4.0*(-1.0)^^(0.5*(3*n + j + 2*m - 2*r).to!double);
					Gamma[row_idx][col_idx] *= sqrt((2.0*n.to!double + 1.0)*(2.0*j.to!double + 1.0));
					Gamma[row_idx][col_idx] /= PI*sqrt(H(m, n)*H(r, j))*(n + j).to!double*(n + j + 2).to!double*((n - j).to!double^^2.0 - 1.0);
				}

				if(((r + m) % 2 == 0) && ((j == (n - 1)) || (j == (n + 1)))) {
					Gamma[row_idx][col_idx] = 1.0;
					Gamma[row_idx][col_idx] /= sqrt(H(m, n)*H(r, j))*sqrt((2.0*n.to!double + 1.0)*(2.0*j.to!double + 1.0));
				}
			},
			// even-even
			(r, j, m, n, row_idx, col_idx) {
				if(r == m) {
					M_c[row_idx][col_idx] = 8.0*(-1.0)^^(0.5*(n + j - 2*m + 2).to!double);
					M_c[row_idx][col_idx] *= sqrt((2.0*n.to!double + 1.0)*(2.0*j.to!double + 1.0));
					M_c[row_idx][col_idx] /= PI*PI*sqrt(H(m, n)*H(m, j))*(n + j).to!double*(n + j + 2).to!double*((n - j).to!double^^2.0 - 1.0);

					D[row_idx][col_idx] = kronecker_delta(j, n)/K(m, n);
				}

				if(((r + m) % 2 != 0) && ((j == (n - 1)) || (j == (n + 1)))) {
					Gamma[row_idx][col_idx] = sgn(r - m);
					Gamma[row_idx][col_idx] /= sqrt(K(m, n)*K(r, j))*sqrt((2.0*n.to!double + 1.0)*(2.0*j.to!double + 1.0));
				}

				if((r + m) % 2 == 0) {
					Gamma[row_idx][col_idx] = 8.0*(-1.0)^^(0.5*(n + j - 2.0*r + 2.0));
					Gamma[row_idx][col_idx] *= sqrt((2*n + 1).to!double*(2*j + 1).to!double);
					Gamma[row_idx][col_idx] /= PI*PI*sqrt(H(m, n)*H(r, j))*(n + j).to!double*(n + j + 2).to!double*((n - j).to!double^^2.0 - 1.0);
				}
			}
		)(Mo, Me);

		//writeln;

		double X = 0;
		iterate_even_odd_sin!(
			(r, j, m, n, row_idx, col_idx) {
				if(r == 0) {
					Meo_s[row_idx][col_idx] = 0;
				} else {
					immutable l = min(r, m);
					immutable exp1 = abs(r - m);
					immutable exp2 = abs(r + m);

					immutable gamma_col_idx = col_idx + odd_states[Mo][0].length;
					immutable gamma_row_idx = row_idx + total_odd_states + even_states[Me][0].length;

					Meo_s[row_idx][col_idx] = (X^^exp1 - ((-1.0)^^l)*X^^exp2)*Gamma[gamma_row_idx][gamma_col_idx];
				}
			}
		)(Mo, Me, 0, 0);

		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, total_even_sin_states, total_odd_sin_states, total_odd_sin_states, 1.0, Meo_s[0].ptr, total_odd_sin_states, A_s_inv[0].ptr, total_odd_sin_states, 0.0, Meo_s_new[0].ptr, total_odd_sin_states);

		iterate_whole_matrix_sin!(
			// odd-odd
			(r, j, m, n, row_idx, col_idx) {

				immutable gamma_col_idx = col_idx + odd_states[Mo][0].length;
				immutable gamma_row_idx = row_idx + odd_states[Mo][0].length;

				if(row_idx == col_idx) {
					M_s[row_idx][col_idx] = 1;
				} else {
					M_s[row_idx][col_idx] = 0;
				}
			},
			// odd-even
			(r, j, m, n, row_idx, col_idx) {
				if(r == 0) {
					M_s[row_idx][col_idx] = 0;
				} else {
					immutable double l = min(r, m);
					immutable double exp1 = abs(r - m);
					immutable double exp2 = abs(r + m);

					immutable gamma_col_idx = col_idx + total_odd_states + even_states[Me][0].length - total_odd_sin_states;
					immutable gamma_row_idx = row_idx + odd_states[Mo][0].length;

					M_s[row_idx][col_idx] = (X^^exp1 - ((-1.0)^^l)*X^^exp2)*Gamma[gamma_row_idx][gamma_col_idx];
				}
			},
			// even-odd
			(r, j, m, n, row_idx, col_idx) {
				M_s[row_idx][col_idx] = Meo_s_new[row_idx - total_odd_sin_states][col_idx];
			},
			// even-even
			(r, j, m, n, row_idx, col_idx) {
				if(r == 0) {
					M_s[row_idx][col_idx] = 0;
				} else {
					immutable l = min(r, m);
					immutable exp1 = abs(r - m);
					immutable exp2 = abs(r + m);

					immutable gamma_col_idx = col_idx + total_odd_states + even_states[Me][0].length - total_odd_sin_states;
					immutable gamma_row_idx = row_idx + total_odd_states + even_states[Me][0].length - total_odd_sin_states;

					M_s[row_idx][col_idx] = (X^^exp1 - ((-1.0)^^l)*X^^exp2)*Gamma[gamma_row_idx][gamma_col_idx];
				}
			}
		)(Mo, Me);

		//Gamma[0][0] 
		// We'll leave it at this size, as this is the same size as the L_c matrix will want.
		ipiv = new int[total_states];

		// We need both M_c and M_c^-1 for timestep updates so copy M_c to M_c_inv and invert.
		foreach(r_idx; 0..total_states) {
			M_c_inv[r_idx][] = M_c[r_idx][];
		}

		foreach(r_idx; 0..total_sin_states) {
			M_s_inv[r_idx][] = M_s[r_idx][];
		}

		iterate_whole_matrix_sin!(
			// odd-odd
			(r, j, m, n, row_idx, col_idx) {
				immutable gamma_col_idx = col_idx + odd_states[Mo][0].length;
				immutable gamma_row_idx = row_idx + odd_states[Mo][0].length;

				D_s[row_idx][col_idx] = D[gamma_row_idx][gamma_col_idx];
			},
			// odd-even
			(r, j, m, n, row_idx, col_idx) {
				immutable gamma_col_idx = col_idx + total_odd_states + even_states[Me][0].length - total_odd_sin_states;
				immutable gamma_row_idx = row_idx + odd_states[Mo][0].length;

				D_s[row_idx][col_idx] = D[gamma_row_idx][gamma_col_idx];
			},
			// even-odd
			(r, j, m, n, row_idx, col_idx) {
				immutable gamma_col_idx = col_idx + odd_states[Mo][0].length;
				immutable gamma_row_idx = row_idx + total_odd_states + even_states[Me][0].length - total_odd_sin_states;

				D_s[row_idx][col_idx] = D[gamma_row_idx][gamma_col_idx];
			},
			// even-even
			(r, j, m, n, row_idx, col_idx) {
				immutable gamma_col_idx = col_idx + total_odd_states + even_states[Me][0].length - total_odd_sin_states;
				immutable gamma_row_idx = row_idx + total_odd_states + even_states[Me][0].length - total_odd_sin_states;

				D_s[row_idx][col_idx] = D[gamma_row_idx][gamma_col_idx];
			}
		)(Mo, Me);

		info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, total_states, total_states, M_c_inv[0].ptr, total_states, ipiv.ptr);
		assert(info == 0, "Failed to invert M_c matrix");
		info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, total_states, M_c_inv[0].ptr, total_states, ipiv.ptr);
		assert(info == 0, "Failed to invert M_c matrix");


		info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, total_sin_states, total_sin_states, M_s_inv[0].ptr, total_sin_states, ipiv.ptr);
		assert(info == 0, "Failed to invert M_s matrix");
		info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, total_sin_states, M_s_inv[0].ptr, total_sin_states, ipiv.ptr);
		assert(info == 0, "Failed to invert M_s matrix");

		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, total_states, total_states, total_states, 1.0, M_c_inv[0].ptr, total_states, D[0].ptr, total_states, 0.0, M_c_inv_D[0].ptr, total_states);

		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, total_sin_states, total_sin_states, total_sin_states, 1.0, M_s_inv[0].ptr, total_sin_states, D_s[0].ptr, total_sin_states, 0.0, M_s_inv_D[0].ptr, total_sin_states);

		iterate_odds!((m, n, idx) {
			K_table[m][n] = K(m, n);
			adjoint_mat[idx] = (-1.0)^^(n.to!double + 1.0);
		})(Mo, 0);

		iterate_evens!((m, n, idx) {
			K_table[m][n] = K(m, n);
			adjoint_mat[idx] = (-1.0)^^(n.to!double + 1.0);
		})(Me, total_odd_states);

		iterate_odds_sin!((m, n, idx) {
			adjoint_mat_sin[idx] = (-1.0)^^(n.to!double + 1.0);
		})(Mo, 0);

		iterate_evens_sin!((m, n, idx) {
			adjoint_mat_sin[idx] = (-1.0)^^(n.to!double + 1.0);
		})(Me, total_odd_sin_states);

		curr_state = 0;
		times[] = 0;
		average_inflow = 0.01;
		chi = 0;
		tau_c[0] = 0.001;
		tau_s[0] = 0.001;

		simple_harmonic_solution(this, 0.0, 0.0);
	}

	@nogc auto find_bracket(double t) {
		
		auto ordered_times = times[get_circular_index(curr_state) + 1..$].chain(times[0..get_circular_index(curr_state) + 1]);
		auto delta_b = time_history - get_circular_index(curr_state);
		ptrdiff_t l = 0;
		ptrdiff_t R = time_history - 1;

		while(l <= R) {
			auto m = (l + R)/2;

			if(m == (time_history - 1)) {
				if(ordered_times[m] != t) {
					return -1;
				} else {
					return m >= delta_b ? m - delta_b + 1 : m - get_circular_index(curr_state) + 1;
				}
			} else if ((ordered_times[m] <= t) && (ordered_times[m + 1] > t)) {
				return m >= delta_b ? m - delta_b + 1 : m - get_circular_index(curr_state) + 1;
			} else if(ordered_times[m] <= t) {
				l = m + 1;
			} else if(ordered_times[m] > t) {
				R = m - 1;
			}
		}

		return -1;
	}

	@nogc auto find_z_bracket(double z) {
		ptrdiff_t l = 0;
		ptrdiff_t R = z_length - 1;

		while(l <= R) {
			auto m = (l + R)/2;

			if(m == (z_length - 1)) {
				return m;
			} else if ((contraction_z_array_alias[m] <= z) && (contraction_z_array_alias[m + 1] > z)) {
				return m;
			} else if(contraction_z_array_alias[m] <= z) {
				l = m + 1;
			} else if(contraction_z_array_alias[m] > z) {
				R = m - 1;
			}
		}

		return 0;
	}

	@nogc void interpolate_to_time(immutable Chunk t, ref Chunk[] time_delay_alpha_buffer, ref Chunk[] time_delay_a_buffer, ref Chunk[] time_delay_beta_buffer, ref Chunk[] time_delay_b_buffer) {

		ptrdiff_t[chunk_size] upper_bounds;
		ptrdiff_t[chunk_size] lower_bounds;

		bool all_nan = t[].map!(a => a.isNaN).fold!((res, a) => res && a)(true);
		if(all_nan) {
			return;
		}

		foreach(idx, ref _t; t[].map!(a => a.isNaN ? 0 : a).enumerate) {
			lower_bounds[idx] = find_bracket(_t);
			upper_bounds[idx] = get_circular_index(lower_bounds[idx] + 1);
		}

		foreach(idx; 0..2*total_states) {
			Chunk t_lower;
			Chunk t_upper;
			Chunk state_u;
			Chunk state_l;
			foreach(c_idx; 0..chunk_size) {
				if(lower_bounds[c_idx] < 0) {
					state_u[c_idx] = 0;
					state_l[c_idx] = 0;
					t_lower[c_idx] = 0;
					t_upper[c_idx] = 1;
				} else {
					t_lower[c_idx] = times[lower_bounds[c_idx]];
					t_upper[c_idx] = times[upper_bounds[c_idx]];

					state_u[c_idx] = state_history[upper_bounds[c_idx]][idx];
					state_l[c_idx] = state_history[lower_bounds[c_idx]][idx];
				}
			}

			immutable Chunk dt = t_upper[] - t_lower[];
			immutable Chunk curr_t = times[get_circular_index(curr_state)];
			immutable Chunk t_delta = curr_t[] - t[];
			immutable Chunk t_delta_abs = abs(t_delta);
			time_delay_alpha_buffer[idx][] = 0;
			
			immutable Chunk time_delta = t[] - t_lower[];
			immutable Chunk state_delta = state_u[] - state_l[];
			time_delay_alpha_buffer[idx][] = state_l[] + time_delta[]*state_delta[]/dt[];

			foreach(c_idx; 0..chunk_size) {
				if(lower_bounds[c_idx] >= 0) {
					if(t_delta_abs[c_idx] <= dt[c_idx]) {
						time_delay_alpha_buffer[idx][c_idx] = state_history[get_circular_index(curr_state)][idx];
					}
				}
			}
		}

		foreach(idx; 0..2*total_sin_states) {
			Chunk t_lower;
			Chunk t_upper;
			Chunk state_u;
			Chunk state_l;
			foreach(c_idx; 0..chunk_size) {
				if(lower_bounds[c_idx] < 0) {
					state_u[c_idx] = 0;
					state_l[c_idx] = 0;
					t_lower[c_idx] = 0;
					t_upper[c_idx] = 1;
				} else {
					t_lower[c_idx] = times[lower_bounds[c_idx]];
					t_upper[c_idx] = times[upper_bounds[c_idx]];

					state_u[c_idx] = state_history[upper_bounds[c_idx]][2*total_states + idx];
					state_l[c_idx] = state_history[lower_bounds[c_idx]][2*total_states + idx];
				}
			}

			immutable Chunk dt = t_upper[] - t_lower[];
			immutable Chunk curr_t = times[get_circular_index(curr_state)];
			immutable Chunk t_delta = curr_t[] - t[];
			immutable Chunk t_delta_abs = abs(t_delta);
			time_delay_beta_buffer[idx][] = 0;
			
			immutable Chunk time_delta = t[] - t_lower[];
			immutable Chunk state_delta = state_u[] - state_l[];
			time_delay_beta_buffer[idx][] = state_l[] + time_delta[]*state_delta[]/dt[];

			foreach(c_idx; 0..chunk_size) {
				if(lower_bounds[c_idx] >= 0) {
					if(t_delta_abs[c_idx] <= dt[c_idx]) {
						time_delay_beta_buffer[idx][c_idx] = state_history[get_circular_index(curr_state)][2*total_states + idx];
					}
				}
			}
		}
		
		foreach(i; 0..total_odd_states) {
			time_delay_a_buffer[i][] = 0;
			time_delay_a_buffer[total_states + i][] = 0;
			foreach(k; 0..total_odd_states) {
				Chunk tmp1 = (time_delay_a_buffer[i][]);
				Chunk tmp2 = A_inv[i][k]*time_delay_alpha_buffer[k][];

				time_delay_a_buffer[i][] = (tmp1[] + tmp2[]);

				tmp1 = (time_delay_a_buffer[total_states + i][]);
				tmp2 = A_inv[i][k]*time_delay_alpha_buffer[total_states + k][];
				time_delay_a_buffer[total_states + i][] = (tmp1[] + tmp2[]);
			}
		}

		foreach(i; 0..total_odd_sin_states) {
			time_delay_b_buffer[i][] = 0;
			time_delay_b_buffer[total_sin_states + i][] = 0;
			foreach(k; 0..total_odd_sin_states) {
				Chunk tmp1 = (time_delay_b_buffer[i][]);
				Chunk tmp2 = A_s_inv[i][k]*time_delay_beta_buffer[k][];

				time_delay_b_buffer[i][] = (tmp1[] + tmp2[]);

				tmp1 = (time_delay_b_buffer[total_sin_states + i][]);
				tmp2 = A_s_inv[i][k]*time_delay_beta_buffer[total_sin_states + k][];
				time_delay_b_buffer[total_sin_states + i][] = (tmp1[] + tmp2[]);
			}
		}

		foreach(i; 0..total_even_states) {
			time_delay_a_buffer[i + total_odd_states][] = time_delay_alpha_buffer[i + total_odd_states][];
			time_delay_a_buffer[i + total_odd_states + total_states][] = time_delay_alpha_buffer[i + total_odd_states + total_states][];
		}

		foreach(i; 0..total_even_sin_states) {
			time_delay_b_buffer[i + total_odd_sin_states][] = time_delay_beta_buffer[i + total_odd_sin_states][];
			time_delay_b_buffer[i + total_odd_sin_states + total_sin_states][] = time_delay_beta_buffer[i + total_odd_sin_states + total_sin_states][];
		}
	}

	@nogc Chunk interpolate_contraction_ratio(immutable Chunk z) {
		Chunk K;

		size_t[chunk_size] lower_bounds;
		size_t[chunk_size] upper_bounds;

		Chunk z_upper;
		Chunk z_lower;
		Chunk K_upper;
		Chunk K_lower;
		foreach(idx, ref _z; z[].map!(a => a.isNaN ? 0 : a).enumerate) {
			lower_bounds[idx] = find_z_bracket(_z);
			if(lower_bounds[idx] == (z_length - 1)) {
				lower_bounds[idx]--;
				upper_bounds[idx] = lower_bounds[idx] + 1;
			} else {
				upper_bounds[idx] = lower_bounds[idx] + 1;
			}

			z_lower[idx] = contraction_z_array_alias[lower_bounds[idx]];
			z_upper[idx] = contraction_z_array_alias[upper_bounds[idx]];

			K_upper[idx] = contraction_array_alias[upper_bounds[idx]];
			K_lower[idx] = contraction_array_alias[lower_bounds[idx]];

			immutable Chunk dz = z_upper[] - z_lower[];
			immutable Chunk z_delta = z[] - z_lower[];
			immutable Chunk K_delta = K_upper[] - K_lower[];
			K[] = (K_lower[] + z_delta[]*K_delta[]/dz[]);
		}

		return K;
	}

	@nogc double wake_skew() {
		return chi;
	}

	@nogc double get_average_inflow() {
		return average_inflow;
	}

	// void update(double C_T, ref RotorInputStateT!AC rotor, ref RotorStateT!AC rotor_state, double advance_ratio, double axial_advance_ratio, ref AircraftStateT!AC ac_state, double dt) {
	// 	update_impl(C_T, rotor, rotor_state, advance_ratio, axial_advance_ratio, ac_state, dt);
	// }

	// void update(double C_T, RotorInputStateT!AC* rotor, RotorStateT!AC* rotor_state, double advance_ratio, double axial_advance_ratio, AircraftStateT!AC* ac_state, double dt) {
	// 	update_impl(C_T, rotor, rotor_state, advance_ratio, axial_advance_ratio, ac_state, dt);
	// }

	void update(ref AircraftInputStateT!AC ac_input , Inflow[] inflows, double freestream_velocity, double advance_ratio, double axial_advance_ratio, double dt) {
		update_impl(ac_input, inflows, freestream_velocity, advance_ratio, axial_advance_ratio, dt);
	}

	void update(AircraftInputStateT!AC* ac_input , Inflow[] inflows, double freestream_velocity, double advance_ratio, double axial_advance_ratio, double dt) {
		update_impl(*ac_input, inflows, freestream_velocity, advance_ratio, axial_advance_ratio, dt);
	}

	package void compute_loading(RS)(auto ref RS rotor_state) {
		tau_c[] = 0;
		tau_s[] = 0;

		immutable omega_sgn = sgn(omega);

		import core.stdc.string : memcpy;
		foreach(b_idx, ref blade_state; rotor_state.blade_states) {

			size_t sin_idx = 0;
			auto _idx = iterate_odds!(
				(m, n, idx) {
					foreach(c_idx, ref chunk; rotor.blades[b_idx].chunks) {
						Chunk atan_num = -omega_sgn*chunk.xi[];
						immutable Chunk psi_r = atan2(atan_num, chunk.r);
						immutable Chunk mpsi = m.to!double*(blade_state.azimuth + psi_r[]);

						Chunk cos_mpsi;
						Chunk sin_mpsi;
						if(m == 0) {
							cos_mpsi[] = 1.0/(2.0*PI);
							sin_mpsi[] = 0;
						} else {
							immutable Chunk[2] sin_cos = sincos(mpsi)[];
							cos_mpsi[] = 1.0/(PI)*sin_cos[1][];
							sin_mpsi[] = 1.0/(PI)*sin_cos[0][];
						}

						Chunk nu = 1.0 - chunk.r[]*chunk.r[];
						nu = sqrt(nu);

						auto blade_frame_forces = Vector!(4, Chunk)(0);
						blade_frame_forces[2][] = blade_state.chunks[c_idx].dC_T[];

						auto global_frame_forces = rotor.blades[b_idx].frame.global_matrix*blade_frame_forces;
						auto rotor_frame_forces = rotor.frame.parent.global_matrix.inverse.get()*global_frame_forces;

						immutable Chunk Pmn = associated_legendre_polynomial_nh(m, n, nu, P_coefficients_nh[idx]);

						//writeln("rotor_frame_forces[2][]: ", rotor_frame_forces[2][]);
						blade_scratch[c_idx][] = rotor_frame_forces[2][]*Pmn[]*cos_mpsi[];

						if(m != 0) {
							blade_scratch_s[c_idx][] = rotor_frame_forces[2][]*Pmn[]*sin_mpsi[];
						}
					}

					tau_c[idx] += integrate_trapaziodal(blade_scratch, rotor.blades[b_idx]);

					if(m != 0) {
						tau_s[sin_idx] += integrate_trapaziodal(blade_scratch_s, rotor.blades[b_idx]);
						sin_idx++;
					}
				}
			)(Mo, 0);

			iterate_evens!(
				(m, n, idx) {
					foreach(c_idx, ref chunk; rotor.blades[b_idx].chunks) {
						Chunk atan_num = -omega_sgn*chunk.xi[];
						immutable Chunk psi_r = atan2(atan_num, chunk.r);
						immutable Chunk mpsi = m.to!double*(blade_state.azimuth + psi_r[]);

						Chunk cos_mpsi;
						Chunk sin_mpsi;
						if(m == 0) {
							cos_mpsi[] = 1.0/(2.0*PI);
							sin_mpsi[] = 0;
						} else {
							immutable Chunk[2] sin_cos = sincos(mpsi)[];
							cos_mpsi[] = 1.0/(PI)*sin_cos[1][];
							sin_mpsi[] = 1.0/(PI)*sin_cos[0][];
						}
						
						Chunk nu = 1.0 - chunk.r[]*chunk.r[];
						nu = sqrt(nu);

						auto blade_frame_forces = Vector!(4, Chunk)(0);
						blade_frame_forces[2][] = blade_state.chunks[c_idx].dC_T[];

						auto global_frame_forces = rotor.blades[b_idx].frame.global_matrix*blade_frame_forces;
						auto rotor_frame_forces = rotor.frame.parent.global_matrix.inverse.get()*global_frame_forces;

						immutable Chunk Pmn = associated_legendre_polynomial(m, n, nu, P_coefficients[idx]);

						blade_scratch[c_idx][] = rotor_frame_forces[2][]*Pmn[]*cos_mpsi[];

						if(m != 0) {
							blade_scratch_s[c_idx][] = rotor_frame_forces[2][]*Pmn[]*sin_mpsi[];
						}
					}

					tau_c[idx] += integrate_trapaziodal(blade_scratch, rotor.blades[b_idx]);

					if(m != 0) {
						tau_s[sin_idx] += integrate_trapaziodal(blade_scratch_s, rotor.blades[b_idx]);
						sin_idx++;
					}
				}
			)(Me, _idx);
		}
	}

	//package void update_impl(RIS, RS, AS)(double C_T, auto ref RIS rotor_input, auto ref RS rotor_state, double _advance_ratio, double _axial_advance_ratio, auto ref AS ac_state, double dt) {
	package void update_impl(ArrayContainer AC = ArrayContainer.None)(ref AircraftInputStateT!AC ac_input , Inflow[] inflows, double freestream_velocity, double _advance_ratio, double _axial_advance_ratio, double dt) {
		omega = rotor_input.angular_velocity;
		
		immutable t_scale = abs(omega);

		advance_ratio = _advance_ratio;
		axial_advance_ratio = -_axial_advance_ratio;

		auto time = dt*curr_state.to!double*t_scale;

		times[get_circular_index(curr_state + 1)] = time;

		compute_loading(rotor_state);

		integrator.step!(system_derivative)(state_history[get_circular_index(curr_state + 1)], state_history[get_circular_index(curr_state)], time, dt*t_scale, this, rotor_input, rotor_state, advance_ratio, axial_advance_ratio);

		alpha = state_history[get_circular_index(curr_state + 1)][0..total_states];
		beta = state_history[get_circular_index(curr_state + 1)][2*total_states..2*total_states + total_sin_states];

		// Get MD state variables
		cblas_dgemv(CblasRowMajor, CblasNoTrans, total_odd_states, total_odd_states, 1.0, A_inv[0].ptr, total_odd_states, alpha.ptr, 1, 0.0, a.ptr, 1);

		cblas_dgemv(CblasRowMajor, CblasNoTrans, total_odd_sin_states, total_odd_sin_states, 1.0, A_s_inv[0].ptr, total_odd_sin_states, beta.ptr, 1, 0.0, b.ptr, 1);

		a[total_odd_states..$] = alpha[total_odd_states..$];
		b[total_odd_sin_states..$] = beta[total_odd_sin_states..$];

		LM_c_scratch[] *= alpha[];

		immutable alpha_0_1 = LM_c_scratch.sum;

		average_inflow = sqrt(3.0)*alpha_0_1;

		foreach(i; 0..total_states) {
			state_history[get_circular_index(curr_state + 1)][total_states + i] = 0;
			foreach(k; 0..total_states) {
				state_history[get_circular_index(curr_state + 1)][total_states + i] += QS_c_mat[i][k]*alpha[k];
			}
		}

		foreach(i; 0..total_sin_states) {
			state_history[get_circular_index(curr_state + 1)][2*total_states + total_sin_states + i] = 0;
			foreach(k; 0..total_sin_states) {
				state_history[get_circular_index(curr_state + 1)][2*total_states + total_sin_states + i] += QS_s_mat[i][k]*beta[k];
			}
		}

		lambda = state_history[get_circular_index(curr_state + 1)][total_states..2*total_states];
		lambda_s = state_history[get_circular_index(curr_state + 1)][2*total_states + total_sin_states..$];

		// Get MD adjoint variables
		cblas_dgemv(CblasRowMajor, CblasNoTrans, total_odd_states, total_odd_states, 1.0, A_inv[0].ptr, total_odd_states, lambda.ptr, 1, 0.0, delta.ptr, 1);

		cblas_dgemv(CblasRowMajor, CblasNoTrans, total_odd_sin_states, total_odd_sin_states, 1.0, A_s_inv[0].ptr, total_odd_sin_states, lambda_s.ptr, 1, 0.0, delta_s.ptr, 1);
		
		delta[total_odd_states..$] = lambda[total_odd_states..$];
		delta_s[total_odd_sin_states..$] = lambda_s[total_odd_sin_states..$];

		//writeln("delta: ", delta, " delta_s: ", delta_s, " lambda: ", lambda, " lambda_s: ", lambda_s, " a: ", a, " b: ", b);

		curr_state++;

		v_0 = compute_inflow_average_at_disk;

		auto v_inf = sqrt(advance_ratio*advance_ratio + axial_advance_ratio*axial_advance_ratio);

		foreach(z_idx, ref z_chunk; contraction_z_array) {
			auto v_z = compute_inflow_average_at_z(z_chunk);
			immutable Chunk v_0_v_z = (v_inf + v_0)/(v_inf + v_z[]);
			contraction_array[z_idx] = sqrt(v_0_v_z);
		}

		global_inverse = local_frame.global_matrix.inverse.get;
		if (advance_ratio > 0) {
			immutable local_freestream = global_inverse*ac_state.freestream;
			immutable normal = Vec4(0, 0, -1, 0);
			immutable projected_freestream = local_freestream - local_freestream.dot(normal)*normal;
			immutable x_axis = Vec4(-1, 0, 0, 0);
			immutable double freestream_rotation = acos(projected_freestream.dot(x_axis)/projected_freestream.magnitude);
			local_frame.rotate(Vec3(0, 0, -1), freestream_rotation);
			local_frame.update(local_frame.parent.global_matrix);
			global_inverse = local_frame.global_matrix.inverse.get;
		}

	}

	double compute_inflow_average_at_disk() {
		Chunk v_z = 0;

		double d_psi = 2.0*PI/n_psi.to!double;
		double d_r = 1.0/n_r.to!double;
		immutable Chunk z = 0;

		foreach(psi_i; 0..n_psi) {
			immutable double psi = 2.0*PI*psi_i.to!double/n_psi.to!double;
			immutable double cos_psi = cos(psi);
			immutable double sin_psi = sin(psi);
			foreach(ref r_idx_chunk; iota(1, n_r+1, 1).chunks(chunk_size)) {
				immutable Chunk r_idx_d = r_idx_chunk.map!(a => a.to!double).staticArray!Chunk;
				immutable Chunk r_chunk = r_idx_d[]/n_r.to!double;
				immutable Chunk x_r = r_chunk[]*cos_psi;
				immutable Chunk y_r = r_chunk[]*sin_psi;
				immutable v = inflow_at_impl(this, x_r, y_r, z);
				v_z[] += v[];
			}
		}

		return v_z.sum*d_r*d_psi*1.0/PI;
	}

	Chunk compute_inflow_average_at_z(immutable Chunk z) {
		Chunk v_z = 0;

		immutable Chunk x_c = -z[]*tan_chi;

		double d_psi = 2.0*PI/n_psi.to!double;
		double d_r = 1.0/n_r.to!double;

		foreach(psi_i; 0..n_psi) {
			immutable double psi = 2.0*PI*psi_i.to!double/n_psi.to!double;
			immutable double cos_psi = cos(psi);
			immutable double sin_psi = sin(psi);
			foreach(r_j; 1..n_r+1) {
				immutable Chunk x_r = (r_j.to!double/n_r.to!double)*cos_psi;
				immutable Chunk y_r = (r_j.to!double/n_r.to!double)*sin_psi;
				immutable Chunk x = x_c[] + x_r[];
				immutable v = inflow_at_impl(this, x, y_r, z);
				v_z[] += v[];
			}
		}

		v_z[] *= d_r*d_psi*1.0/PI;
		return v_z;
	}

	@nogc Chunk compute_contraction_multiplier(immutable Chunk x, immutable Chunk y, immutable Chunk z) {
		immutable K = interpolate_contraction_ratio(z);

		immutable Chunk x_0 = z[]*tan_chi;
		immutable Chunk x_r = x_0[] + x[];

		immutable Chunk r_0_squared = x_r[]*x_r[] + y[]*y[];
		immutable Chunk r_0 = sqrt(r_0_squared);
		
		Chunk k_bar = 0;
		foreach(r_idx, ref _r_0; r_0[]) {
			if(_r_0 <= 1.0) {
				k_bar[r_idx] = K[r_idx];
			} else {
				k_bar[r_idx] = (K[r_idx] + _r_0 - 1.0)/_r_0;
			}
		}

		return k_bar;
	}

	void update_wing_circulation(){
		
	}
	
	Chunk inflow_at(immutable Chunk x, immutable Chunk y, immutable Chunk z, immutable Chunk x_e, double angle_of_attack) {
		if(!contraction_mapping) {
			immutable V = inflow_at_impl(this, x, y, z);
			return V;
		} else {

			immutable k_bar = compute_contraction_multiplier(x, y, z);			

			immutable Chunk x_c = x[]/k_bar[];
			immutable Chunk y_c = y[]/k_bar[];
			immutable V = inflow_at_impl(this, x_c, y_c, z);
			return V;
		}
	}

	@nogc Chunk inflow_at(immutable Chunk r, immutable double cos_azimuth, immutable double sin_azimuth) {
		Chunk i;
		return i;
	}
}
