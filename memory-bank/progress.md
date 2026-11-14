# Progress: OpenCOPTER

## What Works
- ✅ Multirotor aerodynamic simulation (blade element + vortex wake)
- ✅ Huang-Peters dynamic stall inflow model
- ✅ AeroDAS airfoil polar interpolation (XFoil format)
- ✅ Thin airfoil theory model
- ✅ C81 Edgewise dynamic stall model
- ✅ Weissinger L wing lifting surface method
- ✅ Vortex lattice wing method
- ✅ VTK .vtu output (rotors, wings, wakes, wake fields) - needs testing
- ✅ Trim functionality (collective trim to target C_T)
- ✅ BWI noise predictions (broadband + discrete)
- ✅ C API (stable ABI) - 87/125 tests passing
- ✅ C++ API (`include/opencopter.hpp`, `cppbindings.cpp`) - needs test coverage
- ✅ Python bindings via pyd submodule
- ✅ oc_fly Python frontend (JSON config → simulation → post-processing)
- ✅ MPI parallelization across flight conditions
- ✅ AVX/AVX2/AVX512 SIMD optimization via SLEEF
- ✅ HART-II validation cases included
- ✅ Multiple build configurations (debug, release, native, AVX variants)

## What's Left to Build / Known Issues
- ✅ **`oc_blade_geometry_create(NULL airfoil)` fixed** — now defaults to ThinAirfoil theory
- ✅ **`oc_write_rotor_vtu` SEGV fixed (2026-08-01)** — D class FFI casting error resolved. `VtkRotor` is a class (reference type), must cast directly `cast(VtkRotor)` not `cast(VtkRotor*)` + dereference
- ❌ **`oc_wake_create()` causes SEGV** - memory layout mismatch between C/D (crash)
- ⚠️ BladeAirfoil Cl values off by ~2π - may be correct thin airfoil theory, test expectations wrong
- ⚠️ AeroDas GetCd returns 0 for some inputs
- ⚠️ `NullInflow` wake_skew returns NaN instead of 0.0
- ⏳ Windows native support (requires WSL/WSL2 currently)
- ⏳ Additional validation cases beyond HART-II
- ⏳ VTK 9.2+ compatibility (currently locked to 9.1 due to Ubuntu 22.04 repos)
- ⏳ VTK API test coverage needed

## C/C++ Bindings Test Coverage Status
**87 passing / 38 failing / 1 crash** out of 125 tests across 19 test suites

### Passing Categories (✅)
- Frame lifecycle, transforms, hierarchy management
- Aircraft creation/destroy, rotor/wing assignment
- RotorGeometry basic lifecycle
- WingGeometry creation/destroy
- AircraftInputState / RotorInputState all accessors
- AircraftState null-safe operations
- WakeHistory create/destroy/push_back/hybrid
- MemoryBuffer sentinel detection
- Value struct sizes match C expectations

####

