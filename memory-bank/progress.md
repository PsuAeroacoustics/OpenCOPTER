# Progress: OpenCOPTER

## What Works
- **Blade deformation C/C++ API (2026-08-22)**: All 10 blade_inputs fields exposed via C + C++:
  - 4 scalar setters/getters: `set/get_blade_input_{pitch,flapping,flapping_rate,r0}`
  - 6 per-station array setters/getters: `set/get_blade_{flap,lag,twist}_deflection`, `set/get_blade_{flap,lag}_velocity`
  - C header: `include/opencopter.h` (20 new declarations)
  - C++ hpp: `include/opencopter.hpp` (20 new methods on `RotorInputState`)
  - C++ cpp: `source/opencopter/cppbindings.cpp` (20 new implementations)
  - D bindings: `source/opencopter/cbindings.d` (20 new extern(C) functions)
  - 16 new tests in `tests/src/test_blade_input_state.cpp` (BladeInputTestFixture)
  - **175/175 tests pass** (159 original + 16 new)
- **Rotor hub-moment (C_Mx / C_My) tests + ctor bug-fix — P-001-S2 COMPLETE (2026-09-26)**: Added `AircraftState.HubMomentGettersValueAndNullSafe` (C, `test_api_aircraftstate.cpp`: seed `1e9` write-path proof + `isfinite` + `==0.0` + 4 null-safe calls, reusing the `RotorStateCTSetGet` scaffolding) and `CPP_RotorState.HubMomentNullSafe` (C++, `test_cpp_api.cpp`: default-constructed null wrapper → `get_C_Mx()/get_C_My()==0.0`). Both PASS. **Deviation (user-approved):** the test caught a real bug — `RotorStateT` ctors (`aircraft/state.d:195-214`) only zeroed `C_T`, leaving `C_Mx/C_My` NaN pre-step; fixed by adding `C_Q=0; C_Mx=0; C_My=0;` to both ctors. Full suite **180/180 green** (no regressions). Plan **P-001 fully COMPLETE** (S0+S1+S2).
- **Rotor hub-moment (C_Mx / C_My) C++ API — P-001-S1 COMPLETE (2026-09-26)**: Read-only RAII getters `RotorState::get_C_Mx()/get_C_My()` — decls `include/opencopter.hpp:693-694`, one-liner impls `source/opencopter/cppbindings.cpp:761-762` (null-unsafe→0.0, `rst(*this)` helper, mirrors `get_C_T`/`get_C_Q`). Build: `ninja -C build` (recompiled cppbindings.cpp.o) + dub relink clean. `nm -D libopencopter.so` → both C++ symbols + S0 C symbols exported. Unblocks `P-001-S2` (tests).
- **Rotor hub-moment (C_Mx / C_My) C API — P-001-S0 COMPLETE (2026-09-26)**: Read-only getters `oc_rotor_state_get_C_Mx` / `oc_rotor_state_get_C_My` in `include/opencopter.h:363-364` + `source/opencopter/cbindings.d:606-614` (null-guarded, mirroring the `C_T`/`C_Q` pattern). Build clean; both symbols exported from `libopencopter.so`; grep gate: one decl + one def each. Unblocks `P-001-S1` (C++) / `P-001-S2` (tests).

## What's Left
- [x] **Rotor hub-moment (C_Mx / C_My) — P-001 ✅ CLOSED / ARCHIVED (2026-09-26)**: All 3 slices done; 180/180 tests green. Plan at `docs/archive/plans/P-001.md`; all slices at `docs/archive/slices/P-001-S0/S1/S2.md`.
- [ ] **Blade-deformation follow-ups (pre-existing, not from P-001)**:
- [ ] Integration test: 1-step sim with non-zero deflection verifying physics change
- [ ] Python bindings for blade_inputs fields (oc_fly)
- [ ] BLADE_DEFORMATION_PLAN.md final documentation update

## What Worked
# Progress: OpenCOPTER

## What Works
- **Blade deformation C/C++ API (2026-08-22)**: All 10 blade_inputs fields exposed via C + C++:
  - 4 scalar setters/getters: `set/get_blade_input_{pitch,flapping,flapping_rate,r0}`
  - 6 per-station array setters/getters: `set/get_blade_{flap,lag,twist}_deflection`, `set/get_blade_{flap,lag}_velocity`
  - C header: `include/opencopter.h` (20 new declarations)
  - C++ hpp: `include/opencopter.hpp` (20 new methods on `RotorInputState`)
  - C++ cpp: `source/opencopter/cppbindings.cpp` (20 new implementations)
  - D bindings: `source/opencopter/cbindings.d` (20 new extern(C) functions)
  - 16 new tests in `tests/src/test_blade_input_state.cpp` (BladeInputTestFixture)
  - **175/175 tests pass** (159 original + 16 new)
- [ ] Integration test: 1-step sim with non-zero deflection verifying physics change
- [ ] Python bindings for blade_inputs fields (oc_fly)
- [ ] BLADE_DEFORMATION_PLAN.md final documentation update

## What Worked
- ✅ Multirotor aerodynamic simulation (blade element + vortex wake)
- ✅ Huang-Peters dynamic stall inflow model
- ✅ AeroDAS airfoil polar interpolation (XFoil format)
- ✅ Thin airfoil theory model
- ✅ C81 Edgewise dynamic stall model
- ✅ Weissinger L wing lifting surface method
- ✅ Vortex lattice wing method
- ✅ VTK .vtu output (rotors, wings, wakes, wake fields) - needs testing
- ✅ Blade deformation input state (`BladeInputStateT`): per-blade scalars (pitch/flapping/flapping_rate/r_0) + per-station flap/lag/twist deflection + flap/lag velocity arrays; non-breaking `double.infinity` sentinel routing + deflection velocities into `u_p`/`u_t`; `AircraftInputStateT` 4-arg constructor allocates per-station arrays via `num_chunks[]`; C `oc_aircraft_input_state_create_with_chunks` + C++ constructor overload (159/159 tests, 2026-08-22)
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
- ⏳ C/C++/Python bindings for the new `blade_inputs` fields (pitch/flapping/flapping_rate/r_0 + flap/lag/twist deflection) — D-side only so far; legacy arrays remain the settable API
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

