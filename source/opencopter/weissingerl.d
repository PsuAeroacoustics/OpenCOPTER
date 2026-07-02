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

	this(size_t _elements, ref BladeGeometryT!AC blade, double radius, double direction) {
		import std.stdio : writeln;
		
		elements = _elements;
		writeln("elements going in WL model: ", elements);
		immutable chunks = elements/chunk_size;
		writeln("num_chunks going in WL model: ", chunks);
		immutable integration_elements = elements;

		auto influence = allocate_dense(elements, elements);
		//writeln("influence matrix initialized");
		double[][] _influence_inv = allocate_dense(elements, elements);

		//writeln("defined influence_inv");

		immutable m = integration_elements;
		// The collocation stations are psi_v = v*PI/(m + 1) for v = 1 .. m, which
		// is exactly m points. Build them by integer index: a floating-point iota
		// with an upper bound of m*PI/(m + 1) excludes the last station, and its
		// element count is rounding-sensitive across architectures. On AArch64 that
		// left psi_vs/y_array one element short, so y_array[m - 1] read out of
		// bounds (NaN) and poisoned the whole influence matrix.
		immutable psi_vs = iota(1, cast(int)m + 1).map!(v => v*PI/(m + 1.0)).retro.array;
		immutable y_array = psi_vs.map!(psi_mu => cos(psi_mu)).array;
		//writeln("going into nested for loop");
		foreach(ch1; 0..chunks) {
			foreach(c1; 0..chunk_size) {
				immutable v = ch1*chunk_size + c1;

				immutable y = y_array[v];

				immutable true_chord = blade.chunks[ch1].chord[c1]*radius;
				immutable xi_y = direction*blade.chunks[ch1].xi[c1]*radius/true_chord;

				immutable psi_v = psi_vs[v];

				immutable local_aspect = blade.blade_length/(2.0*true_chord);
				//writeln("going into internal nested for loop");
				foreach(ch2; 0..chunks) {
					foreach(c2; 0..chunk_size) {
						immutable n = ch2*chunk_size + c2;

						immutable psi_n = ((elements - n).to!double)*PI/(elements.to!double + 1.0);

						immutable first = 1.0/(elements.to!double + 1.0)*iota(1.0, elements.to!double + 1.0).map!(mu => mu*sin(mu*psi_n)*sin(mu*psi_v)/sin(psi_v)).sum;

						immutable P0 = P(xi_y, -direction*blade.chunks[0].xi[0], y, -1.0, local_aspect)*h_n(elements, psi_n, 0);
						immutable Pend = P(xi_y, -direction*blade.chunks[$-1].xi[$-1], y, 1.0, local_aspect)*h_n(elements, psi_n, PI);

						immutable M = integration_elements.to!double - 1.0;

						immutable Ps = iota(1.0*PI/(M + 1.0), M*PI/(M + 1.0), 1.0*PI/(M + 1.0)).map!((psi_mu) {

							immutable eta_mu = cos(psi_mu);

							immutable eta_idx2 = y_array.countUntil!("a > b")(eta_mu);
							immutable eta_idx1 = eta_idx2 - 1;

							immutable eta_ch_idx1 = eta_idx1/chunk_size;
							immutable eta_c_idx1 = eta_idx1%chunk_size;

							immutable eta_ch_idx2 = eta_idx2/chunk_size;
							immutable eta_c_idx2 = eta_idx2%chunk_size;

							immutable true_chord_eta1 = blade.chunks[eta_ch_idx1].chord[eta_c_idx1]*radius;
							immutable x_eta1 = direction*blade.chunks[eta_ch_idx1].xi[eta_c_idx1]*radius;

							immutable true_chord_eta2 = blade.chunks[eta_ch_idx2].chord[eta_c_idx2]*radius;
							immutable x_eta2 = direction*blade.chunks[eta_ch_idx2].xi[eta_c_idx2]*radius;

							immutable d = y_array[eta_idx2] - y_array[eta_idx1];

							immutable w2 = abs(y_array[eta_idx2] - eta_mu)/d;
							immutable w1 = abs(y_array[eta_idx1] - eta_mu)/d;

							immutable x = w2*x_eta1 + w1*x_eta2;
							immutable true_chord_eta = w2*true_chord_eta1 + w1*true_chord_eta2;

							immutable xi_eta = x/true_chord_eta;

							immutable hn = h_n(elements, psi_n, psi_mu);
							immutable p = P(xi_y, xi_eta, y, eta_mu, local_aspect)*hn;

							return p;
						}).array;

						//writeln("returned P");
						immutable Psum = Ps.sum;
						immutable second = 1.0/(4.0*(M + 1.0))*(
							0.5*(P0 + Pend) + Psum //tehehehehehe
						);

						immutable third = 1.0/(4.0*(M + 1.0))*local_aspect*local_aspect*(

							iota(1.0*PI/(M + 1.0), M*PI/(M + 1.0), 1.0*PI/(M + 1.0)).map!((psi_mu) {

								immutable eta_mu = cos(psi_mu);

								immutable eta_idx2 = y_array.countUntil!("a > b")(eta_mu);
								immutable eta_idx1 = eta_idx2 - 1;
								
								immutable eta_ch_idx1 = eta_idx1/chunk_size;
								immutable eta_c_idx1 = eta_idx1%chunk_size;

								immutable eta_ch_idx2 = eta_idx2/chunk_size;
								immutable eta_c_idx2 = eta_idx2%chunk_size;

								immutable true_chord_eta1 = blade.chunks[eta_ch_idx1].chord[eta_c_idx1]*radius;
								immutable x_eta1 = direction*blade.chunks[eta_ch_idx1].xi[eta_c_idx1]*radius;
								immutable xp_eta1 = blade.chunks[eta_ch_idx1].xi_p[eta_c_idx1]*radius;

								immutable true_chord_eta2 = blade.chunks[eta_ch_idx2].chord[eta_c_idx2]*radius;
								immutable x_eta2 = direction*blade.chunks[eta_ch_idx2].xi[eta_c_idx2]*radius;
								immutable xp_eta2 = blade.chunks[eta_ch_idx2].xi_p[eta_c_idx2]*radius;

								immutable d = y_array[eta_idx2] - y_array[eta_idx1];

								immutable w2 = abs(y_array[eta_idx2] - eta_mu)/d;
								immutable w1 = abs(y_array[eta_idx1] - eta_mu)/d;

								immutable x = w2*x_eta1 + w1*x_eta2;
								immutable true_chord_eta = w2*true_chord_eta1 + w1*true_chord_eta2;
								immutable xi_p_eta = w2*xp_eta1 + w1*xp_eta2;

								immutable xi_eta = x/true_chord_eta;

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
		
		//writeln("calculated influence_inv matrix");
	}

	Chunk compute_bound_circulation_band(BS)(auto ref BS blade_state, size_t chunk_idx, double direction_multiplier, Chunk Cl_alpha, Chunk alpha_zero) {
		Chunk gamma = 0;
		import std.stdio : writeln;


		foreach(c1; 0..chunk_size) {
			immutable r = chunk_idx*chunk_size + c1;
			foreach(ch, ref inf; influence_inv[r]) {
				Chunk tmp = inf[]*sin(blade_state.chunks[ch].aoa)[];
				
				gamma[c1] += tmp.sum;
			}

			gamma[c1] *= sgn(direction_multiplier);
		}
		return gamma;
	}
}
