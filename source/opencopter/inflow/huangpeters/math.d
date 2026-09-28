module opencopter.inflow.huangpeters.math;

import opencopter.math;
import opencopter.memory;

import std.algorithm : min;
import std.conv;
import std.math;
import std.range;

@nogc package double H(long m, long n) {
	return double_factorial(n + m - 1)*double_factorial(n - m - 1)/(double_factorial(n + m)*double_factorial(n - m));
}

@nogc package double K(long m, long n) {
	return ((PI/2.0)^^((-1.0)^^(n.to!double + m.to!double)))*H(m, n);
}

package void zero_matrix(ref double[][] M_c) {
	foreach(ref _M; M_c) {
		_M[] = 0;
	}
}

@nogc package Chunk pow(immutable Chunk x, long power) {
	version(LDC) pragma(inline, true);
	version(GNU) pragma(inline, true);

	Chunk res = 1;
	if(power == 0) {
		res[] = 1;
	} else {
		foreach(idx; 0..power) {
			res[] *= x[];
		}
	}
	return res;
}

@nogc package Chunk associated_legendre_polynomial(bool reduce_order = false)(long m, long n, Chunk x, double[] c) {
	version(LDC) pragma(inline, true);
	version(GNU) pragma(inline, true);

	//x[] = 2.0*x[] - 1.0;

	Chunk p;
	p[] = 0;
	foreach(k; m..n + 1) {
		p[] += c[k - m + 1]*pow(x[], k - m)[];
	}

	immutable Chunk one_m_x2 = 1.0 - x[]*x[];
	Chunk radicand = pow(one_m_x2, m/2)[];
	if((m & 0x1) != 0) {
		// m is odd so we have an integer power plus a square root
		radicand[] *= sqrt(one_m_x2)[];
	}
	p[] *= c[0]*radicand[];
	return p;
}

@nogc package Chunk associated_legendre_polynomial_nh(bool reduce_order = false)(long r, long j, immutable Chunk x, double[] c) {

	version(LDC) pragma(inline, true);
	version(GNU) pragma(inline, true);

	Chunk p = 0;
	Chunk r_bar = 1.0 - x[]*x[];
	r_bar = sqrt(r_bar);
	foreach(qidx, q; iota(r, j, 2).enumerate) {
		p[] += c[q - r + 1]*pow(r_bar[], q)[];
	}
	p[] *= c[0];
	return p;
}

// The m = 0 row of Qmn_bar comes from Q_n's three-term recurrence in n. Q_n is
// the decaying solution of that recurrence, so sweeping it upward amplifies
// rounding by about rho^(2n), rho = eta + sqrt(1 + eta^2), since the growing
// solution P_n takes over. The upward sweep is used only where that stays
// below Q0_upward_growth; further out the row comes from Miller's backward
// sweep, started Q0_miller_depth steps above the top column so that the
// growing solution has decayed below 1e-13 by the time it gets there.
private immutable double Q0_upward_growth = 1.0e3;
private immutable double Q0_miller_log_accuracy = 29.933606208922594; // log(1.0e13)

// Q_n^m(i eta) falls off like eta^-(n + 1). Past this eta every entry is far
// below anything that matters, and evaluating there instead keeps Miller's
// sweep, which grows like eta^n, from overflowing.
private immutable double Q_eta_max = 1.0e15;

/// Smallest eta at which associated_legendre_function takes the m = 0 row of
/// an N-column table from Miller's sweep rather than the upward one.
@nogc package double Q0_miller_eta(long N) {
	return N > 1 ? sinh(log(Q0_upward_growth)/(2.0*(N - 1).to!double)) : double.infinity;
}

@nogc private long Q0_miller_depth(double eta) {
	return cast(long)ceil(Q0_miller_log_accuracy/(2.0*asinh(eta))) + 1;
}

/// K(m, n) for every n > m of an (M + 1) x N table, as associated_legendre_function
/// needs; the other entries are zero.
package double[][] K_table_for(long M, long N) {
	auto K_table = new double[][](M + 1, N);
	foreach(m; 0..M + 1) {
		K_table[m][] = 0;
		foreach(n; m + 1..N) {
			K_table[m][n] = K(m, n);
		}
	}
	return K_table;
}

/// (2n + 1)K(0, n), the coefficients of the m = 0 recurrence, for every n up
/// to the deepest Miller start an N-column table needs.
package double[] Q0_recurrence_coefficients(long N) {
	immutable long top = N - 1 + Q0_miller_depth(Q0_miller_eta(N));
	auto c = new double[top + 1];
	c[0] = 0;
	foreach(n; 1..top + 1) {
		c[n] = (2.0*n.to!double + 1.0)*K(0, n);
	}
	return c;
}

/// Normalised associated Legendre functions of the second kind,
///
///     Qmn_bar[m][n] = Q_n^m(i eta)/Q_n^m(i 0),
///
/// for each lane's eta = x and every column n > m of the table; the diagonal
/// and the entries below it are not needed and are left untouched. K_table[m][n]
/// must hold K(m, n) for every n > m, and Q0_coefficients and Q0_split_eta must
/// come from Q0_recurrence_coefficients and Q0_miller_eta for this table's
/// column count.
@nogc package void associated_legendre_function(Chunk x, ref Chunk[][] Qmn_bar, double[][] K_table, double[] Q0_coefficients, double Q0_split_eta) {

	//version(LDC) pragma(inline, true);
	//version(GNU) pragma(inline, true);

	immutable long N = cast(long)Qmn_bar[0].length;

	// Both m = 0 sweeps run on every lane and each lane keeps the one it takes;
	// NaN lanes take the upward sweep and stay NaN.
	Chunk eta = x;
	double eta_miller_min = double.infinity;
	foreach(i; 0..eta.length) {
		if(eta[i] > Q_eta_max) {
			eta[i] = Q_eta_max;
		}
		if(eta[i] >= Q0_split_eta && eta[i] < eta_miller_min) {
			eta_miller_min = eta[i];
		}
	}

	// std.math.PI is a real, which would take these array expressions out of
	// SIMD registers.
	enum double two_over_pi = 2.0/PI;
	enum double half_pi = PI/2.0;

	// (2/pi) arccot(eta); atan2 keeps full precision where eta is large.
	immutable Chunk one = 1.0;
	immutable Chunk arccot = atan2(one, eta);
	immutable Chunk Q_0 = two_over_pi*arccot[];

	Qmn_bar[0][0][] = Q_0[];
	Qmn_bar[0][1][] = 1.0 - half_pi*eta[]*Qmn_bar[0][0][];
	foreach(n; 1..N - 1) {
		immutable double c = Q0_coefficients[n];
		Qmn_bar[0][n + 1][] = Qmn_bar[0][n - 1][] - c*eta[]*Qmn_bar[0][n][];
	}

	if(eta_miller_min < double.infinity) {
		// Miller's algorithm: sweep down from zero above the top column, deep
		// enough for the smallest eta taking it, then scale to the exact Q_0.
		immutable long top = min(N - 1 + Q0_miller_depth(eta_miller_min), cast(long)Q0_coefficients.length - 1);
		Chunk q_above = 0.0;
		Chunk q = 1.0;
		long n = top;
		// Nothing is kept above column N - 1, so rescale every few steps to
		// keep lanes with large eta in range on a deep sweep.
		for(; n > N; n--) {
			immutable double c = Q0_coefficients[n];
			immutable Chunk q_below = q_above[] + c*eta[]*q[];
			q_above = q;
			q = q_below;
			if((top - n)%8 == 7) {
				q_above[] /= q[];
				q[] = 1.0;
			}
		}
		q_above[] /= q[];
		q[] = 1.0;
		for(; n > 0; n--) {
			immutable double c = Q0_coefficients[n];
			immutable Chunk q_below = q_above[] + c*eta[]*q[];
			q_above = q;
			q = q_below;
			foreach(i; 0..eta.length) {
				if(eta[i] >= Q0_split_eta) {
					Qmn_bar[0][n - 1][i] = q[i];
				}
			}
		}
		foreach(i; 0..eta.length) {
			if(eta[i] >= Q0_split_eta) {
				immutable double Q_0_over_q = Q_0[i]/q[i];
				foreach(k; 0..N) {
					Qmn_bar[0][k][i] *= Q_0_over_q;
				}
			}
		}
	}

	// m > 0 from
	//
	//     Q_n^(m+1) = (Q_(n-1)^m - (n - m) K(m, n) eta Q_n^m)/sqrt(1 + eta^2),
	//
	// which follows from sqrt(z^2 - 1) Q_n^(m+1) = (n - m) z Q_n^m - (n + m) Q_(n-1)^m
	// at z = i eta. It stays accurate for any eta once row 0 is, and entry n
	// of row m + 1 needs only entries n - 1 and n of row m, so starting each
	// row at n = m + 2 fills every n > m.
	immutable Chunk xs = eta[]*eta[];
	immutable Chunk one_plus_xs = 1.0 + xs[];
	immutable Chunk sqrt_opxs = sqrt(one_plus_xs);
	immutable Chunk one_over_sqrt_opxs = 1.0/sqrt_opxs[];

	foreach(m; 0..cast(long)Qmn_bar.length - 1) {
		foreach(n; m + 2..N) {
			immutable Knmm = (n.to!double - m.to!double)*K_table[m][n];
			immutable Chunk Qtmp = Knmm*eta[]*Qmn_bar[m][n][];
			immutable Chunk Qmn_delta = Qmn_bar[m][n - 1][] - Qtmp[];
			Qmn_bar[m + 1][n][] = one_over_sqrt_opxs[]*Qmn_delta[];
		}
	}
}

unittest {
	// Reference values from mpmath, legenq(n, m, i eta, type=3)/legenq(n, m, i 0+, type=3),
	// on the largest table the state lists produce. The etas cover both m = 0
	// sweeps and the far field where the old upward sweep had lost every digit.
	immutable long M = 12;
	immutable long N = 15;
	auto K_table = K_table_for(M, N);
	auto Qmn_bar = new Chunk[][](M + 1, N);
	Chunk eta = [0.0, 0.05, 0.3, 1.0, 5.0, 100.0, 1.0e4, 1.0e6];
	associated_legendre_function(eta, Qmn_bar, K_table, Q0_recurrence_coefficients(N), Q0_miller_eta(N));

	static struct Expected { long m; long n; double[8] Q; }
	immutable Expected[] expected = [
		Expected(0, 1, [1.0, 9.2395810344635232e-1, 6.1619814030489117e-1, 2.1460183660255169e-1, 1.3022200750596208e-2, 3.3331333476179366e-5, 3.3333333133333333e-9, 3.3333333333313331e-13]),
		Expected(0, 14, [1.0, 4.8396092865040313e-1, 1.3416763036140309e-2, 2.3545917636927096e-6, 1.2040976125344985e-15, 4.2784109662991539e-35, 4.2800672168893146e-65, 4.2800673825527704e-95]),
		Expected(1, 2, [1.0, 8.8829232576250872e-1, 4.8608147563998405e-1, 1.0168585115698146e-1, 1.5426160367269961e-3, 1.9998143022604175e-7, 1.9999999814285715e-13, 1.9999999999981429e-19]),
		Expected(3, 8, [1.0, 6.7085206326264313e-1, 9.1991032104022138e-2, 6.6945927432965878e-4, 2.0910136243132628e-9, 4.5237082284181305e-21, 4.5248867599428459e-39, 4.5248868778162658e-57]),
		Expected(6, 11, [1.0, 6.1118953795305087e-1, 5.22960963092177e-2, 1.0264890663841485e-4, 5.7570490025986144e-12, 1.6338935776118483e-27, 1.6345210225563926e-51, 1.6345210853157242e-75]),
		Expected(12, 13, [1.0, 7.2663999297435045e-1, 1.3087789627923493e-1, 5.2685392858246838e-4, 4.7760281146471706e-12, 3.7014439684737021e-30, 3.7037034776500718e-58, 3.7037037036810983e-86]),
		Expected(12, 14, [1.0, 6.617933957283616e-1, 7.7602418780400806e-2, 1.3793906634864241e-4, 2.6963413084380408e-13, 1.0485783354418517e-32, 1.0492278819425596e-62, 1.0492279469204994e-92]),
	];
	foreach(e; expected) {
		foreach(i; 0..eta.length) {
			assert(isClose(Qmn_bar[e.m][e.n][i], e.Q[i], 1.0e-12, 0.0));
		}
	}
}

@nogc Chunk sign(Chunk x) {
	Chunk res;
	foreach(idx, ref _x; x) {
		if(_x > 0.0) {
			res[idx] = 1.0;
		} else {
			res[idx] = -1.0;
		}
	}
	return res;
}

@nogc double sign(double x) {
	if(x > 0.0) {
		return 1.0;
	} else {
		return -1.0;
	}
}
