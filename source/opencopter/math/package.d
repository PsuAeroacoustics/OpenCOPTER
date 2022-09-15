module opencopter.math;

public import opencopter.math.algebraic;
public import opencopter.math.complex;
public import opencopter.math.trigonometry;
public import opencopter.math.integration;

import opencopter.aircraft;
import opencopter.math.sleef;

import numd.linearalgebra.matrix : Matrix, Vector;

import core.simd;

import std.algorithm : map;
import std.traits : isStaticArray, Unqual, isFloatingPoint, ForeachType, isSIMDVector;

alias Vec3 = Vector!(3, double);
alias Mat4 = Matrix!(4, 4, double);

void print_matlab(double[][] mat) {
	import std.stdio : writeln, write;
	write("[");
	foreach(ref r; mat[0..$-1]) {
		writeln(r, ";");
	}
	write(mat[$-1]);
	writeln("]");
}

double[] sinspace(double start, double end, size_t elements) {
	
	import std.algorithm : map;
	import std.array : array;
	import std.conv : to;
	import std.math : sin, PI;
	import std.range : iota, zip;

	double R = end - start;

	immutable d_theta = (0.5*PI)/(elements + 1).to!double;
	auto thetas = iota(0, 0.5*PI, d_theta);
	auto r = thetas.map!(a => R*sin(a));

	return zip(r[0..$-1], r[1..$]).map!(a => (a[0] + 0.5*(a[1] - a[0])).to!double).array;
}

double[] cosspace(double start, double end, size_t elements) {
	
	import std.algorithm : map;
	import std.array : array;
	import std.conv : to;
	import std.math : cos, PI;
	import std.range : iota, zip;

	double R = end - start;

	immutable d_theta = (0.5*PI)/(elements + 1).to!double;
	auto thetas = iota(0, 0.5*PI, d_theta);
	auto r = thetas.map!(a => R*cos(a));

	return zip(r[0..$-1], r[1..$]).map!(a => (a[0] + 0.5*(a[1] - a[0])).to!double).array;
}

version(LDC) {
	pragma(LDC_intrinsic, "llvm.sqrt.f#")
    @nogc T llvm_rsqrt(T)(T val)
        if (__traits(isFloating, T));

}

@nogc T rsqrt(T)(auto ref T vector) if(isStaticArray!T || isSIMDVector!T) {
	static import std.math;
	Unqual!T result;
	version(LDC) {
		import ldc.intrinsics : llvm_sqrt;
		double8 vec = cast(double8)vector[];
		result = cast(T)llvm_rsqrt(vec);
	} else {
		foreach(idx; 0..T.length) {
			result[idx] = 1.0/std.math.sqrt(vector[idx]);
		}
	}

	return result;
}

@nogc T sqrt(T)(auto ref T vector) if(isStaticArray!T || isSIMDVector!T) {
	static import std.math;
	Unqual!T result;

	static if(use_sleef) {
		version(AVX_512F) {
			version(LDC) {
				static if(T.length == 16) {
					import ldc.intrinsics : llvm_sqrt;
					double8 vec1 = cast(double8)vector[0..8];
					double8 vec2 = cast(double8)vector[8..$];
					result[0..8] = cast(double[8])llvm_sqrt(vec1);
					result[8..$] = cast(double[8])llvm_sqrt(vec2);
				} else {
					import ldc.intrinsics : llvm_sqrt;
					double8 vec = cast(double8)vector[];
					result = cast(T)llvm_sqrt(vec);
				}
				
			} else {
				foreach(idx; 0..T.length) {
					result[idx] = std.math.sqrt(vector[idx]);
				}
			}
		} else {
			version(D_AVX2) {
				pragma(msg, "avx2 sqrt");
				static if(T.length == 16) {
					result[0..4] = cast(double[4])Sleef_sqrtd4_avx2(cast(double4)vector[0..4]);
					result[4..8] = cast(double[4])Sleef_sqrtd4_avx2(cast(double4)vector[4..8]);
					result[8..12] = cast(double[4])Sleef_sqrtd4_avx2(cast(double4)vector[8..12]);
					result[12..$] = cast(double[4])Sleef_sqrtd4_avx2(cast(double4)vector[12..$]);
				} else static if(T.length == 8) {
					result[0..4] = cast(double[4])Sleef_sqrtd4_avx2(cast(double4)vector[0..4]);
					result[4..$] = cast(double[4])Sleef_sqrtd4_avx2(cast(double4)vector[4..$]);
				} else static if(T.length == 4) {
					result[] = cast(T)Sleef_sqrtd4_avx2(cast(double4)vector);
				}
			} else {
				version(D_AVX) {
					static if(T.length == 8) {
						result[0..2] = cast(double[2])Sleef_sqrtd2_sse4(cast(double2)vector[0..2]);
						result[2..4] = cast(double[2])Sleef_sqrtd2_sse4(cast(double2)vector[2..4]);
						result[4..6] = cast(double[2])Sleef_sqrtd2_sse4(cast(double2)vector[4..6]);
						result[6..$] = cast(double[2])Sleef_sqrtd2_sse4(cast(double2)vector[6..$]);

					} else static if(T.length == 4) {
						result[0..2] = cast(double[2])Sleef_sqrtd2_sse4(cast(double2)vector[0..2]);
						result[2..$] = cast(double[2])Sleef_sqrtd2_sse4(cast(double2)vector[2..$]);
					} else static if(T.length == 2) {
						result[] = cast(T)Sleef_sqrtd2_sse4(cast(double2)vector);
					}

				} else {
					version(X86_64) {
						static if(T.length == 8) {
							result[0..2] = cast(double[2])Sleef_sqrtd2_sse2(cast(double2)vector[0..2]);
							result[2..4] = cast(double[2])Sleef_sqrtd2_sse2(cast(double2)vector[2..4]);
							result[4..6] = cast(double[2])Sleef_sqrtd2_sse2(cast(double2)vector[4..6]);
							result[6..$] = cast(double[2])Sleef_sqrtd2_sse2(cast(double2)vector[6..$]);

						} else static if(T.length == 4) {
							result[0..2] = cast(double[2])Sleef_sqrtd2_sse2(cast(double2)vector[0..2]);
							result[2..$] = cast(double[2])Sleef_sqrtd2_sse2(cast(double2)vector[2..$]);
						} else static if(T.length == 2) {
							result[] = cast(T)Sleef_sqrtd2_sse2(cast(double2)vector);
						}
					} else {
						foreach(idx; 0..T.length) {
							result[idx] = std.math.sqrt(vector[idx]);
						}
					}
					
				}
			}
		}
	} else {
		foreach(idx; 0..T.length) {
			result[idx] = std.math.sqrt(vector[idx]);
		}
	}
	return result;
}

@nogc T erf(T)(auto ref T vector) if(isStaticArray!T) {
	static import std.mathspecial;
	Unqual!T result;

	import std.stdio : writeln;

	version(LDC) {
		version(D_AVX2) {
			pragma(msg, "avx2 erf");
			static if(T.length == 8) {
				result[0..4] = cast(double[4])Sleef_erfd4_u10avx2(cast(double4)vector[0..4]);
				result[4..$] = cast(double[4])Sleef_erfd4_u10avx2(cast(double4)vector[4..$]);
			} else static if(T.length == 4) {
				result[] = cast(T)Sleef_erfd4_u10avx2(cast(double4)vector);
			}
		} else {
			version(D_AVX) {
				static if(T.length == 8) {
					result[0..2] = cast(double[2])Sleef_erfd2_u10sse4(cast(double2)vector[0..2]);
					result[2..4] = cast(double[2])Sleef_erfd2_u10sse4(cast(double2)vector[2..4]);
					result[4..6] = cast(double[2])Sleef_erfd2_u10sse4(cast(double2)vector[4..6]);
					result[6..$] = cast(double[2])Sleef_erfd2_u10sse4(cast(double2)vector[6..$]);

				} else static if(T.length == 4) {
					result[0..2] = cast(double[2])Sleef_erfd2_u10sse4(cast(double2)vector[0..2]);
					result[2..$] = cast(double[2])Sleef_erfd2_u10sse4(cast(double2)vector[2..$]);
				} else static if(T.length == 2) {
					result[] = cast(T)Sleef_erfd2_u10sse4(cast(double2)vector);
				}

			} else {
				version(X86_64) {
					static if(T.length == 8) {
						result[0..2] = cast(double[2])Sleef_erfd2_u10sse2(cast(double2)vector[0..2]);
						result[2..4] = cast(double[2])Sleef_erfd2_u10sse2(cast(double2)vector[2..4]);
						result[4..6] = cast(double[2])Sleef_erfd2_u10sse2(cast(double2)vector[4..6]);
						result[6..$] = cast(double[2])Sleef_erfd2_u10sse2(cast(double2)vector[6..$]);

					} else static if(T.length == 4) {
						result[0..2] = cast(double[2])Sleef_erfd2_u10sse2(cast(double2)vector[0..2]);
						result[2..$] = cast(double[2])Sleef_erfd2_u10sse2(cast(double2)vector[2..$]);
					} else static if(T.length == 2) {
						result[] = cast(T)Sleef_erfd2_u10sse2(cast(double2)vector);
					}
				} else {
					foreach(idx; 0..T.length) {
						result[idx] = std.mathspecial.erf(vector[idx]);
					}
				}
				
			}
		}
	} else {
		foreach(idx; 0..T.length) {
			result[idx] = std.mathspecial.erf(vector[idx]);
		}
	}
	return result;
}


@nogc T abs(T)(auto ref T vector) if(isStaticArray!T || isSIMDVector!T) {
	static import std.math;
	Unqual!T result;
	static foreach(idx; 0..T.length) {
		result[idx] = std.math.abs(vector[idx]);
	}

	return result;
}

@nogc T sgn(T)(auto ref T vector) if(isStaticArray!T) {
	static import std.math;
	Unqual!T result;
	static foreach(idx; 0..T.length) {
		result[idx] = std.math.sgn(vector[idx]);
	}

	return result;
}

@nogc ForeachType!T sum(T)(T v) if(isStaticArray!T) {
	import std.array : staticArray;
	alias F = ForeachType!T;
	static if(T.length == 2) {
		return v[0] + v[1];
	} else {
		F[T.length/2] r = v[0..T.length/2] + v[T.length/2..T.length];
		return sum(r);
	}
}

/+import opencopter.memory : Chunk;
import std.math : floor;
@nogc Chunk sum(Chunk[] data) {
	if(data.length <= 2) {
		Chunk s = 0;
		foreach(ref d; data) {
			s[] += d[];
		}
		return s;
	} else {
		uint m = cast(uint)floor(cast(double)data.length/2.0);
		immutable Chunk s0 = sum(data[0..m]);
		immutable Chunk s1 = sum(data[m..$]);
		immutable Chunk s = s0[] + s1[];
		return s;
	}
}+/
