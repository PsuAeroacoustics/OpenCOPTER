module opencopter.weissingerl;

import opencopter.aircraft;
import opencopter.config;
import opencopter.math;
import opencopter.math.blas;
import opencopter.math.lapacke;
import opencopter.memory;
import opencopter.io;

import numd.utility;

import std.algorithm;
import std.conv;
import std.math;
import std.range;

// Nitya: CIRCULATION MODEL

private auto P(double xi_y, double xi_eta, double y, double eta, double local_aspect) {
	immutable global_divisor = 1.0/(y - eta);
	immutable numerator = xi_y - xi_eta + 0.5;
	immutable denominator = sqrt(
		(xi_y - xi_eta + 0.5)^^2.0
		+ local_aspect*local_aspect*(y - eta)^^2.0
	);
	immutable arg = numerator/denominator - 1;
	if(abs(arg) < 1.0e-14) {
		return 0;
	} else {
		return global_divisor*arg;
	}
}

private auto R(double xi_y, double xi_eta, double xi_p_eta, double y, double eta, double local_aspect) {
	immutable numerator = xi_y - xi_eta + 0.5 + xi_p_eta*(eta - y);
	immutable denominator = (
		(xi_y - xi_eta + 0.5)^^2.0
		+ local_aspect*local_aspect*(y - eta)^^2.0
	)^^1.5;

	return numerator/denominator;
}

private double h_n(size_t m, double psi_n, double psi_v) {
	import std.algorithm : map, sum;
	import std.range : iota;
	return 2.0/(m.to!double + 1.0)*iota(1.0, m.to!double + 1).map!(mu => mu*sin(mu*psi_n)*cos(mu*psi_v)).sum;
}

private double f_n(size_t m, double psi_n, double psi_v) {
	import std.algorithm : map, sum;
	import std.range : iota;
	return 2.0/(m.to!double + 1.0)*iota(1.0, m.to!double + 1).map!(mu => sin(mu*psi_n)*sin(mu*psi_v)).sum;
}

struct WeissingerL(ArrayContainer AC) {

	private Chunk[][] influence_inv;
	private size_t elements;

	this(size_t _elements, ref BladeGeometryT!AC blade, double radius) {
		import std.stdio : writeln;
		
		elements = _elements;
		immutable chunks = elements/chunk_size;

		immutable integration_elements = elements;

		auto influence = allocate_dense(elements, elements);
		double[][] _influence_inv = allocate_dense(elements, elements);

		auto y_array = generate_radius_points(elements);

		foreach(ch1; 0..chunks) {
			foreach(c1; 0..chunk_size) {
				immutable v = ch1*chunk_size + c1;
				immutable y = 2.0*y_array[v] - 1.0;

				immutable true_chord = blade.chunks[ch1].chord[c1]*radius;
				immutable xi_y = blade.chunks[ch1].xi[c1]*radius/true_chord;

				immutable psi_v = ((elements - v).to!double)*PI/(elements.to!double + 1.0);

				immutable local_aspect = blade.blade_length/(2.0*true_chord);

				foreach(ch2; 0..chunks) {
					foreach(c2; 0..chunk_size) {
						immutable n = ch2*chunk_size + c2;

						immutable psi_n = ((elements - n).to!double)*PI/(elements.to!double + 1.0);

						immutable first = 1.0/(elements.to!double + 1.0)*iota(1.0, elements.to!double + 1.0).map!(mu => mu*sin(mu*psi_n)*sin(mu*psi_v)/sin(psi_v)).sum;

						immutable second = 1.0/(4.0*(integration_elements.to!double + 1.0))*(
							0.5*(P(xi_y, blade.chunks[0].xi[0], y, -1.0, local_aspect)*h_n(elements, psi_n, 0) + P(xi_y, blade.chunks[$-1].xi[$-1], y, 1.0, local_aspect)*h_n(elements, psi_n, PI)) +
							iota(1.0, integration_elements.to!double + 1.0).enumerate.map!((e) {
								immutable mu = e[1];
								immutable eta_idx = e[0];
								immutable psi_mu = mu*PI/(integration_elements.to!double + 1.0);
								immutable eta_mu = cos(psi_mu);
								immutable eta_ch_idx = eta_idx/chunk_size;
								immutable eta_c_idx = eta_idx%chunk_size;

								immutable true_chord_eta = blade.chunks[eta_ch_idx].chord[eta_c_idx]*radius;
								immutable xi_eta = blade.chunks[eta_ch_idx].xi[eta_c_idx]*radius/true_chord_eta;

								immutable p = P(xi_y, xi_eta, y, eta_mu, local_aspect)*h_n(elements, psi_n, psi_mu);

								return p;
							}).array.sum
						);

						immutable third = 1.0/(4.0*(integration_elements.to!double + 1.0))*local_aspect*local_aspect*(
							iota(1.0, integration_elements.to!double + 1.0).enumerate.map!((e) {
								immutable mu = e[1];
								immutable eta_idx = e[0];
								immutable psi_mu = mu*PI/(integration_elements.to!double + 1.0);
								immutable eta_mu = cos(psi_mu);

								immutable eta_ch_idx = eta_idx/chunk_size;
								immutable eta_c_idx = eta_idx%chunk_size;

								immutable true_chord_eta = blade.chunks[eta_ch_idx].chord[eta_c_idx]*radius;

								immutable xi_eta = blade.chunks[eta_ch_idx].xi[eta_c_idx]*radius/true_chord_eta;
								immutable xi_p_eta = blade.chunks[eta_ch_idx].xi_p[eta_c_idx];///true_chord_eta;

								return R(xi_y, xi_eta, xi_p_eta, y, eta_mu, local_aspect)*f_n(elements, psi_n, psi_mu)*sin(psi_mu);
							}).array.sum
						);

						// The 0.5 is to compensate for the fact that we are
						// actually integrating over half the length as the
						// original formulation.
						influence[v][n] = 0.5*(first - second + third);
					}
				}
			}
		}

		foreach(r_idx; 0..elements) {
			_influence_inv[r_idx][] = influence[r_idx][];
		}

		debug writeln;
		debug _influence_inv.print_matlab;
		openblas_set_num_threads(1);

		int info = 0;
		auto ipiv = new int[elements];
		info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, elements.to!int, elements.to!int, _influence_inv[0].ptr, elements.to!int, ipiv.ptr);
		assert(info == 0, "Failed to invert influence matrix");
		info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, elements.to!int, _influence_inv[0].ptr, elements.to!int, ipiv.ptr);

		assert(info == 0, "Failed to invert influence matrix");

		influence_inv = allocate_dense_chunk_aliased(elements, elements);

		foreach(r; 0..elements) {
			foreach(ch; 0..chunks) {
				influence_inv[r][ch][] = _influence_inv[r][ch*chunk_size..ch*chunk_size+chunk_size];
			}
		}
	}

	Chunk compute_bound_circulation_band(BS)(auto ref BS blade_state, size_t chunk_idx, double direction_multiplier, Chunk Cl_alpha, Chunk alpha_zero) {
		Chunk gamma = 0;
		import std.stdio : writeln;

		//writeln("Cl_alpha: ", Cl_alpha, ", alpha_zero: ", alpha_zero);
		
		foreach(c1; 0..chunk_size) {
			immutable r = chunk_idx*chunk_size + c1;
			foreach(ch, ref inf; influence_inv[r]) {
				//Chunk aoa_eff = Cl_alpha[]/(2.0*PI)*(blade_state.chunks[ch].aoa[] - alpha_zero[]);
				//Chunk tmp = inf[]*sin(aoa_eff)[];

				//writeln("aoa_eff: ", aoa_eff, ", blade_state.chunks[ch].aoa: ", blade_state.chunks[ch].aoa);
				Chunk tmp = inf[]*sin(blade_state.chunks[ch].aoa)[];
				
				gamma[c1] += tmp.sum;
			}

			gamma[c1] *= -sgn(direction_multiplier);
		}
		return gamma;
	}
}
