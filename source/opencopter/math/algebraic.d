module opencopter.math.algebraic;

import std.conv;

@nogc double factorial(long n) {
	double r = 1;
	foreach(i; 1..n + 1) {
		r *= i.to!double;
	}
	return r;
}

unittest {
	assert(0.factorial == 1);
	assert(1.factorial == 1);
	assert(2.factorial == 2);
	assert(3.factorial == 6);
	assert(4.factorial == 24);
	assert(5.factorial == 120);
	assert(6.factorial == 720);
	assert(7.factorial == 5_040);
	assert(8.factorial == 40_320);
	assert(9.factorial == 362_880);
	assert(10.factorial == 3_628_800);
}

@nogc double falling_factorial(double x, long n) {
	double r = 1;
	foreach(k; 0..n) {
		r *= (x - k.to!double);
	}
	return r;
}

unittest {
	assert(2.falling_factorial(0) == 1);
	assert(2.falling_factorial(1) == 2);
	assert(2.falling_factorial(2) == 2);
	assert(2.falling_factorial(3) == 0);
	assert(2.falling_factorial(4) == 0);
	
	assert(5.falling_factorial(0) == 1);
	assert(5.falling_factorial(1) == 5);
	assert(5.falling_factorial(2) == 20);
}

@nogc double double_factorial(long n) {
	
	if(n < 0) {
		return (n + 2).double_factorial/(n + 2).to!double;
	} else if(n == 0) {
		return 1.0;
	} else if(n % 2 == 0) {
		double fac = 2.0;
		foreach(k; 2..n/2 + 1) {
			fac *= 2.0*k.to!double;
		}
		return fac;
	} else if(n % 2 != 0) {
		double fac = 1;
		foreach(k; 1..(n + 1)/2 + 1) {
			fac *= 2.0*k.to!double - 1.0;
		}
		return fac;
	}
	assert(false);
}

unittest {
	assert((-5).double_factorial == 1./3.);
	assert((-3).double_factorial == -1);
	assert((-2).double_factorial == double.infinity);
	assert((-1).double_factorial == 1);
	assert(0.double_factorial == 1);
	assert(1.double_factorial == 1);
	assert(2.double_factorial == 2);
	assert(3.double_factorial == 3);
	assert(4.double_factorial == 8);
	assert(5.double_factorial == 15);
	assert(6.double_factorial == 48);
	assert(7.double_factorial == 105);
	assert(8.double_factorial == 384);
	assert(9.double_factorial == 945);
}

@nogc double kronecker_delta(long i, long j) {
	if(i == j) {
		return 1;
	}
	return 0;
}

unittest {
	assert(kronecker_delta(0, 0) == 1);
	assert(kronecker_delta(0, 1) == 0);
	assert(kronecker_delta(1, 1) == 1);
	assert(kronecker_delta(1, 2) == 0);
	assert(kronecker_delta(2, 2) == 1);
}

/+@nogc T pow(T)(auto ref T vector, immutable ForeachType!T power) if(isStaticArray!T) {

	static import std.math;
	Unqual!T result;

	import std.stdio : writeln;

	version(LDC) {
		version(D_AVX2) {
			pragma(msg, "avx2 cos");
			static if(T.length == 8) {
				result[0..4] = cast(double[4])Sleef_powd4_u10avx2(cast(double4)vector[0..4]);
				result[4..$] = cast(double[4])Sleef_powd4_u10avx2(cast(double4)vector[4..$]);
			} else static if(T.length == 4) {
				result[] = cast(T)Sleef_powd4_u10avx2(cast(double4)vector);
			}
		} else {
			version(D_AVX) {
				static if(T.length == 8) {
					result[0..2] = cast(double[2])Sleef_powd2_u10sse4(cast(double2)vector[0..2]);
					result[2..4] = cast(double[2])Sleef_powd2_u10sse4(cast(double2)vector[2..4]);
					result[4..6] = cast(double[2])Sleef_powd2_u10sse4(cast(double2)vector[4..6]);
					result[6..$] = cast(double[2])Sleef_powd2_u10sse4(cast(double2)vector[6..$]);

				} else static if(T.length == 4) {
					result[0..2] = cast(double[2])Sleef_cosd2_u10sse4(cast(double2)vector[0..2]);
					result[2..$] = cast(double[2])Sleef_cosd2_u10sse4(cast(double2)vector[2..$]);
				} else static if(T.length == 2) {
					result[] = cast(T)Sleef_powd2_u10sse4(cast(double2)vector);
				}

			} else {
				version(X86_64) {
					static if(T.length == 8) {
						result[0..2] = cast(double[2])Sleef_cosd2_u10sse2(cast(double2)vector[0..2]);
						result[2..4] = cast(double[2])Sleef_cosd2_u10sse2(cast(double2)vector[2..4]);
						result[4..6] = cast(double[2])Sleef_cosd2_u10sse2(cast(double2)vector[4..6]);
						result[6..$] = cast(double[2])Sleef_cosd2_u10sse2(cast(double2)vector[6..$]);

					} else static if(T.length == 4) {
						result[0..2] = cast(double[2])Sleef_cosd2_u10sse2(cast(double2)vector[0..2]);
						result[2..$] = cast(double[2])Sleef_cosd2_u10sse2(cast(double2)vector[2..$]);
					} else static if(T.length == 2) {
						result[] = cast(T)Sleef_cosd2_u10sse2(cast(double2)vector);
					}
				} else {
					foreach(idx; 0..T.length) {
						result[idx] = std.math.cos(vector[idx]);
					}
				}
				
			}
		}
	} else {
		foreach(idx; 0..T.length) {
			result[idx] = std.math.cos(vector[idx]);
		}
	}
	return result;
}
+/