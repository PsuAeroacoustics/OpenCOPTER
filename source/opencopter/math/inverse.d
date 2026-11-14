module opencopter.math.inverse;

import opencopter.math.blas : openblas_set_num_threads;
import opencopter.math.lapacke;

import std.conv : to;
import std.exception : enforce;
import std.format : format;
import std.math : abs, fmax;

/++
 +	Largest residual tolerated by invert_checked, as max |a_inv*a - I|.
 +	A correct inverse of the influence matrices we build lands near 1e-13;
 +	the wrong inverses seen on CI land near 1e2.
 +/
enum double inverse_residual_tolerance = 1.0e-6;

/++
 +	Inverts `a` into `a_inv` with LAPACK getrf/getri, then checks the result.
 +
 +	LAPACK has been seen on CI to report success (info == 0) while returning an
 +	inverse that is wrong in every entry, so `info` alone can't be trusted. The
 +	check multiplies the result back against `a` with plain D loops, not BLAS,
 +	so a broken BLAS kernel can't hide its own error.
 +
 +	On failure it throws with two extra facts that tell the likely causes apart:
 +	whether `a` changed while LAPACK ran (something wrote over our memory), and
 +	whether inverting a fresh copy of `a` reproduces the same wrong answer (the
 +	library computes it wrong, deterministically, on this machine).
 +
 +	`a_inv` must be an n x n matrix whose rows are one contiguous row-major
 +	block, as allocate_dense makes them. Its contents are overwritten.
 +/
void invert_checked(double[][] a_inv, const double[][] a, string name) {
	immutable n = a.length;
	enforce(a_inv.length == n, name~": inverse has "~a_inv.length.to!string~" rows, matrix has "~n.to!string);
	foreach(r; 0..n) {
		enforce(a[r].length == n && a_inv[r].length == n, name~": matrix is not square");
		enforce(a_inv[r].ptr == a_inv[0].ptr + r*n, name~": inverse rows are not one contiguous block");
	}

	auto snapshot = new double[n*n];
	foreach(r; 0..n) {
		snapshot[r*n..(r + 1)*n] = a[r][];
		a_inv[r][] = a[r][];
	}

	lu_invert(a_inv[0].ptr, n, name);

	immutable residual = inverse_residual(a_inv, a);
	if(residual <= inverse_residual_tolerance) return;

	size_t changed = 0;
	foreach(r; 0..n) {
		foreach(c; 0..n) {
			// Compare bits, so a NaN that stayed NaN doesn't count as a change.
			if(a[r][c] !is snapshot[r*n + c]) changed++;
		}
	}

	auto retry_block = snapshot.dup;
	lu_invert(retry_block.ptr, n, name);
	auto retry = new double[][](n);
	auto original = new double[][](n);
	foreach(r; 0..n) {
		retry[r] = retry_block[r*n..(r + 1)*n];
		original[r] = snapshot[r*n..(r + 1)*n];
	}
	immutable retry_residual = inverse_residual(retry, original);
	bool retry_identical = true;
	foreach(r; 0..n) {
		if(retry[r] != a_inv[r]) {
			retry_identical = false;
			break;
		}
	}

	throw new Exception(format!(
		"%s: LAPACK returned a wrong %sx%s inverse with info == 0 (max |inv*A - I| = %.3g, tolerance %.0g). "~
		"Input entries changed while LAPACK ran: %s. Inverting a fresh copy gave max |inv*A - I| = %.3g, %s the first result."
	)(name, n, n, residual, inverse_residual_tolerance, changed, retry_residual,
		retry_identical ? "bit-identical to" : "different from"));
}

private void lu_invert(double* a, size_t n, string name) {
	openblas_set_num_threads(1);
	auto ipiv = new int[n];
	int info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n.to!int, n.to!int, a, n.to!int, ipiv.ptr);
	enforce(info == 0, format!"%s: LU factorisation failed (getrf info = %s)"(name, info));
	info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, n.to!int, a, n.to!int, ipiv.ptr);
	enforce(info == 0, format!"%s: inversion failed (getri info = %s)"(name, info));
}

private double inverse_residual(const double[][] a_inv, const double[][] a) {
	immutable n = a.length;
	auto row = new double[n];
	double worst = 0;
	foreach(r; 0..n) {
		row[] = 0;
		foreach(k; 0..n) {
			row[] += a_inv[r][k]*a[k][];
		}
		row[r] -= 1.0;
		foreach(v; row) {
			// fmax drops NaN, so test for it explicitly.
			if(v != v) return double.nan;
			worst = fmax(worst, abs(v));
		}
	}
	return worst;
}

unittest {
	import opencopter.memory : allocate_dense;

	// A diagonally dominant matrix inverts cleanly and passes the check.
	immutable n = 24;
	auto a = allocate_dense(n, n);
	foreach(r; 0..n) {
		foreach(c; 0..n) {
			a[r][c] = r == c ? n + 1.0 : 1.0/(1.0 + r + 2.0*c);
		}
	}
	auto a_inv = allocate_dense(n, n);
	invert_checked(a_inv, a, "test matrix");
	assert(inverse_residual(a_inv, a) < 1.0e-12);

	// A wrong inverse is caught by the residual.
	a_inv[3][5] += 1.0e-3;
	assert(inverse_residual(a_inv, a) > inverse_residual_tolerance);
}
