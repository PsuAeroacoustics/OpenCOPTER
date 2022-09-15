module opencopter.math.trigonometry;

import std.traits;

import opencopter.math.sleef;

import core.simd;

@nogc T cos(T)(auto ref T vector) if(isStaticArray!T) {

	static import std.math;
	Unqual!T result;

	import std.stdio : writeln;

	//version(LDC) {
	static if(use_sleef) {
		version(D_AVX2) {
			pragma(msg, "avx2 cos");
			static if(T.length == 16) {
				result[0..4] = cast(double[4])Sleef_cosd4_u10avx2(cast(double4)vector[0..4]);
				result[4..8] = cast(double[4])Sleef_cosd4_u10avx2(cast(double4)vector[4..8]);
				result[8..12] = cast(double[4])Sleef_cosd4_u10avx2(cast(double4)vector[8..12]);
				result[12..$] = cast(double[4])Sleef_cosd4_u10avx2(cast(double4)vector[12..$]);
			} else static if(T.length == 8) {
				result[0..4] = cast(double[4])Sleef_cosd4_u10avx2(cast(double4)vector[0..4]);
				result[4..$] = cast(double[4])Sleef_cosd4_u10avx2(cast(double4)vector[4..$]);
			} else static if(T.length == 4) {
				result[] = cast(T)Sleef_cosd4_u10avx2(cast(double4)vector);
			}
		} else {
			version(D_AVX) {
				pragma(msg, "avx cos");
				static if(T.length == 16) {
					result[0..2] = cast(double[2])Sleef_cosd2_u10sse4(cast(double2)vector[0..2]);
					result[2..4] = cast(double[2])Sleef_cosd2_u10sse4(cast(double2)vector[2..4]);
					result[4..6] = cast(double[2])Sleef_cosd2_u10sse4(cast(double2)vector[4..6]);
					result[6..8] = cast(double[2])Sleef_cosd2_u10sse4(cast(double2)vector[6..8]);
					result[8..10] = cast(double[2])Sleef_cosd2_u10sse4(cast(double2)vector[8..10]);
					result[10..12] = cast(double[2])Sleef_cosd2_u10sse4(cast(double2)vector[10..12]);
					result[12..14] = cast(double[2])Sleef_cosd2_u10sse4(cast(double2)vector[12..14]);
					result[14..$] = cast(double[2])Sleef_cosd2_u10sse4(cast(double2)vector[14..$]);
				} else static if(T.length == 8) {
					result[0..2] = cast(double[2])Sleef_cosd2_u10sse4(cast(double2)vector[0..2]);
					result[2..4] = cast(double[2])Sleef_cosd2_u10sse4(cast(double2)vector[2..4]);
					result[4..6] = cast(double[2])Sleef_cosd2_u10sse4(cast(double2)vector[4..6]);
					result[6..$] = cast(double[2])Sleef_cosd2_u10sse4(cast(double2)vector[6..$]);

				} else static if(T.length == 4) {
					result[0..2] = cast(double[2])Sleef_cosd2_u10sse4(cast(double2)vector[0..2]);
					result[2..$] = cast(double[2])Sleef_cosd2_u10sse4(cast(double2)vector[2..$]);
				} else static if(T.length == 2) {
					result[] = cast(T)Sleef_cosd2_u10sse4(cast(double2)vector);
				}

			} else {
				version(X86_64) {
					pragma(msg, "sse cos");
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

@nogc T sin(T)(auto ref T vector) if(isStaticArray!T) {
	import std.math : sin;
	Unqual!T result;
	foreach(idx, ref v; vector) {
		result[idx] = sin(v);
	}
	return result;
}

@nogc T tan(T)(auto ref T vector) if(isStaticArray!T) {
	static import std.math;
	Unqual!T result;
	foreach(idx; 0..T.length) {
		result[idx] = std.math.tan(vector[idx]);
	}
	return result;
}

@nogc T atan(T)(auto ref T vector) if(isStaticArray!T) {
	static import std.math;
	Unqual!T result;

	import std.stdio : writeln;

	//version(LDC) {
	static if(use_sleef) {
		version(D_AVX2) {
			pragma(msg, "avx2 atan");
			static if(T.length == 16) {
				result[0..4] = cast(double[4])Sleef_atand4_u10avx2(cast(double4)vector[0..4]);
				result[4..8] = cast(double[4])Sleef_atand4_u10avx2(cast(double4)vector[4..8]);
				result[8..12] = cast(double[4])Sleef_atand4_u10avx2(cast(double4)vector[8..12]);
				result[12..$] = cast(double[4])Sleef_atand4_u10avx2(cast(double4)vector[12..$]);
			} else static if(T.length == 8) {
				result[0..4] = cast(double[4])Sleef_atand4_u10avx2(cast(double4)vector[0..4]);
				result[4..$] = cast(double[4])Sleef_atand4_u10avx2(cast(double4)vector[4..$]);
			} else static if(T.length == 4) {
				result[] = cast(T)Sleef_atand4_u10avx2(cast(double4)vector);
			}
		} else {
			version(D_AVX) {
				static if(T.length == 8) {
					result[0..2] = cast(double[2])Sleef_atand2_u10sse4(cast(double2)vector[0..2]);
					result[2..4] = cast(double[2])Sleef_atand2_u10sse4(cast(double2)vector[2..4]);
					result[4..6] = cast(double[2])Sleef_atand2_u10sse4(cast(double2)vector[4..6]);
					result[6..$] = cast(double[2])Sleef_atand2_u10sse4(cast(double2)vector[6..$]);

				} else static if(T.length == 4) {
					result[0..2] = cast(double[2])Sleef_atand2_u10sse4(cast(double2)vector[0..2]);
					result[2..$] = cast(double[2])Sleef_atand2_u10sse4(cast(double2)vector[2..$]);
				} else static if(T.length == 2) {
					result[] = cast(T)Sleef_atand2_u10sse4(cast(double2)vector);
				}

			} else {
				version(X86_64) {
					static if(T.length == 8) {
						result[0..2] = cast(double[2])Sleef_atand2_u10sse2(cast(double2)vector[0..2]);
						result[2..4] = cast(double[2])Sleef_atand2_u10sse2(cast(double2)vector[2..4]);
						result[4..6] = cast(double[2])Sleef_atand2_u10sse2(cast(double2)vector[4..6]);
						result[6..$] = cast(double[2])Sleef_atand2_u10sse2(cast(double2)vector[6..$]);

					} else static if(T.length == 4) {
						result[0..2] = cast(double[2])Sleef_atand2_u10sse2(cast(double2)vector[0..2]);
						result[2..$] = cast(double[2])Sleef_atand2_u10sse2(cast(double2)vector[2..$]);
					} else static if(T.length == 2) {
						result[] = cast(T)Sleef_atand2_u10sse2(cast(double2)vector);
					}
				} else {
					foreach(idx; 0..T.length) {
						result[idx] = std.math.atan(vector[idx]);
					}
				}
				
			}
		}
	} else {
		foreach(idx; 0..T.length) {
			result[idx] = std.math.atan(vector[idx]);
		}
	}
	return result;
}

@nogc T atan2(T)(auto ref T vector1, auto ref T vector2) if(isStaticArray!T) {
	static import std.math;
	Unqual!T result;

	import std.stdio : writeln;

	//version(LDC) {
	static if(use_sleef) {
		version(D_AVX2) {
			pragma(msg, "avx2 atan2");
			static if(T.length == 16) {
				result[0..4] = cast(double[4])Sleef_atan2d4_u10avx2(cast(double4)vector1[0..4], cast(double4)vector2[0..4]);
				result[4..8] = cast(double[4])Sleef_atan2d4_u10avx2(cast(double4)vector1[4..8], cast(double4)vector2[4..8]);
				result[8..12] = cast(double[4])Sleef_atan2d4_u10avx2(cast(double4)vector1[8..12], cast(double4)vector2[8..12]);
				result[12..$] = cast(double[4])Sleef_atan2d4_u10avx2(cast(double4)vector1[12..$], cast(double4)vector2[12..$]);
			} else static if(T.length == 8) {
				result[0..4] = cast(double[4])Sleef_atan2d4_u10avx2(cast(double4)vector1[0..4], cast(double4)vector2[0..4]);
				result[4..$] = cast(double[4])Sleef_atan2d4_u10avx2(cast(double4)vector1[4..$], cast(double4)vector2[4..$]);
			} else static if(T.length == 4) {
				result[] = cast(T)Sleef_atan2d4_u10avx2(cast(double4)vector1, cast(double4)vector2);
			}
		} else {
			version(D_AVX) {
				static if(T.length == 8) {
					result[0..2] = cast(double[2])Sleef_atan2d2_u10sse4(cast(double2)vector1[0..2], cast(double2)vector2[0..2]);
					result[2..4] = cast(double[2])Sleef_atan2d2_u10sse4(cast(double2)vector1[2..4], cast(double2)vector2[2..4]);
					result[4..6] = cast(double[2])Sleef_atan2d2_u10sse4(cast(double2)vector1[4..6], cast(double2)vector2[4..6]);
					result[6..$] = cast(double[2])Sleef_atan2d2_u10sse4(cast(double2)vector1[6..$], cast(double2)vector2[6..$]);

				} else static if(T.length == 4) {
					result[0..2] = cast(double[2])Sleef_atan2d2_u10sse4(cast(double2)vector1[0..2], cast(double2)vector2[0..2]);
					result[2..$] = cast(double[2])Sleef_atan2d2_u10sse4(cast(double2)vector1[2..$], cast(double2)vector2[2..$]);
				} else static if(T.length == 2) {
					result[] = cast(T)Sleef_atan2d2_u10sse4(cast(double2)vector1, cast(double2)vector2);
				}

			} else {
				version(X86_64) {
					static if(T.length == 8) {
						result[0..2] = cast(double[2])Sleef_atan2d2_u10sse2(cast(double2)vector1[0..2], cast(double2)vector2[0..2]);
						result[2..4] = cast(double[2])Sleef_atan2d2_u10sse2(cast(double2)vector1[2..4], cast(double2)vector2[2..4]);
						result[4..6] = cast(double[2])Sleef_atan2d2_u10sse2(cast(double2)vector1[4..6], cast(double2)vector2[4..6]);
						result[6..$] = cast(double[2])Sleef_atan2d2_u10sse2(cast(double2)vector1[6..$], cast(double2)vector2[6..$]);

					} else static if(T.length == 4) {
						result[0..2] = cast(double[2])Sleef_atan2d2_u10sse2(cast(double2)vector1[0..2], cast(double2)vector2[0..2]);
						result[2..$] = cast(double[2])Sleef_atan2d2_u10sse2(cast(double2)vector1[2..$], cast(double2)vector2[2..$]);
					} else static if(T.length == 2) {
						result[] = cast(T)Sleef_atan2d2_u10sse2(cast(double2)vector1, cast(double2)vector2);
					}
				} else {
					foreach(idx; 0..T.length) {
						result[idx] = std.math.atan2(vector1[idx], vector2[idx]);
					}
				}
				
			}
		}
	} else {
		foreach(idx; 0..T.length) {
			result[idx] = std.math.atan2(vector1[idx], vector2[idx]);
		}
	}
	return result;
}
