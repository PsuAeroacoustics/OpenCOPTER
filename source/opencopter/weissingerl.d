module opencopter.weissingerl;

import opencopter.aircraft;
import opencopter.config;
import opencopter.math;
import opencopter.math.blas;
import opencopter.math.lapacke;
import opencopter.memory;

import numd.utility;

import std.algorithm;
import std.conv;
import std.math;
import std.range;

private auto P(double sweep, double y, double eta, double aspect) {
	immutable tan_sweep = tan(sweep);
	immutable global_divisor = 1.0/(y - eta);
	immutable numerator = y*tan_sweep - eta*tan_sweep + 0.5;
	immutable denominator = sqrt(
		(y*tan_sweep - eta*tan_sweep + 0.5)^^2.0
		+ 0.25*aspect*aspect*(y - eta)^^2.0
	);
	immutable arg = numerator/denominator - 1;
	if(abs(arg) < 1.0e-14) {
		return 0;
	} else {
		return global_divisor*arg;
	}
}

private auto R(double sweep, double y, double eta, double aspect) {
	immutable tan_sweep = tan(sweep);
	immutable numerator = y*tan_sweep - eta*tan_sweep + 0.5 + tan_sweep*(y - eta);
	immutable denominator = (
		(y*tan_sweep - eta*tan_sweep + 0.5)^^2.0
		+ 0.25*aspect*aspect*(y - eta)^^2.0
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

	this(size_t _elements, ref BladeGeometryT!AC blade) {
		import std.stdio : writeln;
		
		elements = _elements;
		immutable chunks = elements/chunk_size;

		immutable integration_elements = elements;

		auto influence = allocate_dense(elements, elements);
		double[][] _influence_inv = allocate_dense(elements, elements);

		foreach(ch1; 0..chunks) {
			foreach(c1; 0..chunk_size) {
				immutable v = ch1*chunk_size + c1;
				immutable y = 2.0*blade.chunks[ch1].r[c1] - 1.0;

				immutable psi_v = ((elements - v).to!double)*PI/(elements.to!double + 1.0);
				foreach(ch2; 0..chunks) {
					foreach(c2; 0..chunk_size) {
						immutable n = ch2*chunk_size + c2;
						
						immutable eta = 2.0*blade.chunks[ch2].r[c2] - 1.0;
						
						immutable sweep = blade.chunks[ch1].sweep[c1];
						immutable aspect = 1.0/blade.chunks[ch1].chord[c1];

						immutable psi_n = ((elements - n).to!double)*PI/(elements.to!double + 1.0);

						immutable first = 1.0/(elements.to!double + 1.0)*iota(1.0, elements.to!double + 1.0).map!(mu => mu*sin(mu*psi_n)*sin(mu*psi_v)/sin(psi_v)).sum;

						immutable second = 1.0/(4.0*(integration_elements.to!double + 1.0))*(
							0.5*(P(sweep, y, -1.0, aspect)*h_n(elements, psi_n, 0) + P(sweep, y, 1.0, aspect)*h_n(elements, psi_n, PI)) +
							iota(1.0, integration_elements.to!double + 1.0).map!((mu) {
								immutable psi_mu = mu*PI/(integration_elements.to!double + 1.0);
								immutable eta_mu = cos(psi_mu);

								return P(sweep, y, eta_mu, aspect)*h_n(elements, psi_n, psi_mu);
							}).array.sum
						);

						immutable third = 1.0/(16.0*(integration_elements.to!double + 1.0))*aspect*aspect*(
							iota(1.0, integration_elements.to!double + 1.0).map!((mu) {
								immutable psi_mu = mu*PI/(integration_elements.to!double + 1.0);
								immutable eta_mu = cos(psi_mu);

								return R(sweep, y, eta_mu, aspect)*f_n(elements, psi_n, psi_mu)*sin(psi_mu);
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

	void compute_bound_circulation(BS)(auto ref BS blade_state) {
		foreach(ch1; 0..influence_inv[0].length) {
			foreach(c1; 0..chunk_size) {
				immutable r = ch1*chunk_size + c1;
				double gamma = 0;
				foreach(ch, ref inf; influence_inv[r]) {
					Chunk tmp = inf[]*sin(blade_state.chunks[ch].aoa)[];
					gamma += tmp.sum;
				}

				blade_state.chunks[ch1].d_gamma[c1] = blade_state.chunks[ch1].gamma[c1] - gamma;
				blade_state.chunks[ch1].gamma[c1] = gamma;
			}
		}
	}

	void compute_bound_circulation_band(BS)(auto ref BS blade_state, immutable Chunk u, size_t chunk_idx) {
		foreach(c1; 0..chunk_size) {
			immutable r = chunk_idx*chunk_size + c1;
			double gamma = 0;
			foreach(ch, ref inf; influence_inv[r]) {
				Chunk tmp = inf[]*sin(blade_state.chunks[ch].aoa)[];
				gamma += tmp.sum;
			}

			gamma *= /+5*+/u[c1];
			blade_state.chunks[chunk_idx].d_gamma[c1] = blade_state.chunks[chunk_idx].gamma[c1] - gamma;
			//blade_state.chunks[chunk_idx].d_gamma[c1] = gamma - blade_state.chunks[chunk_idx].gamma[c1];
			blade_state.chunks[chunk_idx].gamma[c1] = gamma;
		}
	}
}
