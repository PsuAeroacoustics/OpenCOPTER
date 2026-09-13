module opencopter.inflow.huangpeters.math;

import opencopter.math;
import opencopter.memory;

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

@nogc package void associated_legendre_function(Chunk x, ref Chunk[][] Qmn_bar, double[][] K_table) {

	//version(LDC) pragma(inline, true);
	//version(GNU) pragma(inline, true);

	immutable Chunk xs = x[]*x[];

	immutable Chunk one_plus_xs = 1.0 + xs[];
	immutable Chunk sqrt_opxs = sqrt(one_plus_xs);
	immutable Chunk one_over_sqrt_opxs = 1.0/sqrt_opxs[];

	// m = 0 row: Q_n(x), x = eta >= 0.
	//
	// Q_n(x) for x > 0 is the *minimal* (decaying) solution of the three-term
	// recurrence, so an upward sweep is numerically unstable (roundoff is
	// amplified by the dominant, growing solution).  Instead we use Miller's
	// backward recurrence: seed the top with arbitrary values, sweep downward
	// (which suppresses the dominant solution), then rescale to the exact base
	// value Q_0(x).
	//
	// The base value uses the cancellation-free identity
	//     PI/2 - atan(x) = atan(1/x)   (x >= 0)
	// so no precision is lost for large x.
	immutable long N0 = Qmn_bar[0].length - 1;
	if(N0 >= 2) {
		Qmn_bar[0][N0][] = 1.0;
		Qmn_bar[0][N0 - 1][] = 1.0;
		// Backward form of  Q[n+1] = Q[n-1] - (2n+1) K[0][n] x Q[n].
		for(long n = N0 - 1; n >= 1; n--) {
			immutable Kc1 = (2.0*n.to!double + 1.0)*K_table[0][n];
			Qmn_bar[0][n - 1][] = Qmn_bar[0][n + 1][] + Kc1*x[]*Qmn_bar[0][n][];
		}
		// Rescale the whole row so that Q[0] matches the exact base value.
		immutable Chunk inv_x = 1.0/x[];
		immutable Chunk Qtrue0 = (2.0/PI)*atan(inv_x)[];
		immutable Chunk scale = Qtrue0[]/Qmn_bar[0][0][];
		foreach(n; 0..Qmn_bar[0].length) {
			Qmn_bar[0][n][] = scale[]*Qmn_bar[0][n][];
		}
	} else {
		// Tiny-N edge case: use the exact closed forms directly.
		immutable Chunk inv_x0 = 1.0/x[];
		Qmn_bar[0][0][] = (2.0/PI)*atan(inv_x0)[];
		if(N0 >= 1) {
			Qmn_bar[0][1][] = 1.0 - x[]*PI/2.0*Qmn_bar[0][0][];
		}
	}

	// m >= 1 rows.
	//
	// The original applied the (1+x^2)^{-1/2} factor per step in a *separate*
	// pass after the whole m-recurrence, leaving intermediate values unbounded.
	// We instead fold the normalisation into the recurrence so each row is
	// computed already in its final (bounded) scale.
	//
	// The original's deferred normalisation yields, for row index r >= 1,
	//     Q[r] = f^(r-1) U[r],   f = 1/sqrt(1+x^2),
	// which is algebraically equivalent to the per-step recurrence
	//     Q[1]   = Q[0][n-1] - n K[0][n] x Q[0][n]           (no factor)
	//     Q[r+1] = f (Q[r][n-1] - (n-r) K[r][n] x Q[r][n])   (r >= 1)
	// i.e. the factor is applied when stepping from row m to m+1 for m >= 1.
	foreach(m; 0..Qmn_bar.length - 1) {
		foreach(n; 1..Qmn_bar[m + 1].length) {
			immutable Knmm = (n.to!double - m.to!double)*K_table[m][n];
			immutable Chunk Qtmp = Knmm*x[]*Qmn_bar[m][n][];
			immutable Chunk Qmn_delta = Qmn_bar[m][n - 1][] - Qtmp[];
			if(m >= 1) {
				Qmn_bar[m + 1][n][] = one_over_sqrt_opxs[]*Qmn_delta[];
			} else {
				Qmn_bar[m + 1][n][] = Qmn_delta[];
			}
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
