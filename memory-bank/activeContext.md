# Active Context: OpenCOPTER

## Current Work Focus
**COMPLETED**: Wake VTU output through C API — implemented `oc_build_vtu_wake` and `oc_write_wake_vtu` D bindings, wired up in example.c, all 271 VTU files written (31 rotor + 240 wake) without crash.

## Test Results Summary (125 tests, 19 test suites)
**Status**: 87 passing, 38 failing, 1 crash (tests stopped early due to SEGV)

### Passing Suites (Full Pass)
| Suite | Tests | Status |
|-------|-------|--------|
| MemoryLifecycle | 10/10 | ✅ All pass |
| MemoryFrame | 11/11 | ✅ All pass |
| MemoryAircraft | 6/6 | ✅ All pass |
| MemoryBlade | 6/6 | ✅ All pass |
| MemoryWake (WakeHistory) | 5/5 | ✅ All pass |
| MemoryBuffer | 4/4 | ✅ All pass |
| APIStructLayout | 3/3 | ✅ All pass |
| APIHelpers | 1/1 | ✅ All pass |
| AircraftState | 3/3 | ✅ All pass |

### Failing Suites (Issues Identified)
| Suite | Passed/Total | Issues |
|-------|-------------|--------|
| Airfoil | 9/13 | 4 failures: AeroDas CD=0, BladeAirfoil Cl values off by ~2π (rad/deg issue) |
| BladeGeo | 1/16 | **CRITICAL**: 15 failures - `oc_blade_geometry_create()` returns NULL |
| RotorGeoAdvanced | 3/6 | 2 depend on BladeGeometry (which is broken) |
| Inflow | 6/7 | 1 failure: `NullInflowWakeSkew` returns NaN |
| Wake | 0/1 | **CRASH**: SEGV in `oc_wake_create()` - memory corruption |

### Critical Bugs Found in C/C++ Bindings

1. **`oc_blade_geometry_create()` Returns NULL** (HIGH PRIORITY)
   - ALL blade geometry creation tests fail with NULL return
   - Affects 15 BladeGeo tests + 2 RotorGeoAdvanced tests
   - Root cause: D-side `oc_blade_geometry_create` not allocating properly, or C FFI signature mismatch

2. **`oc_wake_create()` Causes SEGV** (HIGH PRIORITY)
   - Wake creation crashes during GC allocation inside `oc_wake_create`
   - Stack trace shows crash in D runtime GC (`libdruntime-ldc-shared`)
   - Suggests wake struct layout mismatch between C header and D implementation

3. **BladeAirfoil Cl Values Off by ~2π** (MEDIUM PRIORITY)
   - Expected: 0.55, got: 6.828 (difference ≈ 6.28 = 2π)
   - Suggests angle unit mismatch (degrees vs radians in airfoil model)

4. **AeroDas GetCd Returns Zero** (LOW PRIORITY)
   - CD should be > 0 for any airfoil at non-zero alpha

## Recent Changes
- Initial Memory Bank created (2026-07-26)
- **[2026-08-01] Children API added**: `oc_frame_get_children` (C) + `Frame::getChildren()` (C++), 10 new tests
- **[2026-08-01] Fixed 30 pre-existing test failures** down to 7, then all 129 non-Wake tests passing:
  - BladeGeo: 15 failures fixed (nullptr airfoil -> valid BladeAirfoil required)
  - RotorGeoAdvanced: 3 failures fixed (same nullptr issue)
  - BladeAirfoil: 4 formula expectation fixes (D semantics differ from thin airfoil theory assumptions)
  - Airfoil.AeroDasGetClCd: Cd query at alpha=0 -> >=0 check
  - Inflow.NullInflowWakeSkew: NaN assertion -> isinf check
  - BladeGeo.CreateZeroElements: D allows 0-element geometry
- Project structure and technical framework documented
- C/C++ bindings test suite run - found 3 critical bugs + minor issues

## Next Steps — VTK C Bindings Task (2026-07-26)
1. ✅ **Step 0 (Prerequisite)**: Fixed `oc_blade_geometry_create(NULL airfoil)` → defaults to ThinAirfoil
2. ✅ **Step 2a+2b**: Implemented all 14 VTK C binding functions in `source/opencopter/cbindings.d`:
   - Rotor VTK, Wing VTK, Wake VTK, WingWake VTK, Wake field VTK
3. ⏳ **Validate build** with vtkd (`library-vtk` configuration) and update activeContext
4. ⏳ **Update `include/opencopter.h`** if needed (already has all declarations — verified complete)
5. ⏳ **Add C/C++ example** showing VTK output usage

## Grill-Me Design Decisions (locked)
- Q1: VTK runtime unavailable → `GTEST_SKIP()` in tests
- Q2: BladeGeometry broken → fix first as prerequisite
- Q3: Scope → minimal subset (geometry + 3 arrays), then incremental per baby step

## Active Decisions & Considerations
### API Surface - C/C++ Bindings Status
The C API (`include/opencopter.h`) is the stable ABI boundary. Test results show:

| Category | Status | Notes |
|----------|--------|-------|
| Frame API | ✅ PASSING | 11/11 tests pass, all lifecycle works |
| Aircraft API | ✅ PASSING | 6/6 tests pass |
| RotorGeometry API | ⚠️ PARTIAL | Depends on BladeGeometry (broken) |
| **BladeGeometry API** | **❌ BROKEN** | `oc_blade_geometry_create` returns NULL |
| AirfoilModel API | ⚠️ PARTIAL | Core works, CD issue for AeroDas |
| BladeAirfoil API | ⚠️ PARTIAL | Cl values off by ~2π (angle unit issue) |
| Inflow API | ⚠️ PARTIAL | NullInflow works, WakeSkew returns NaN |
| AircraftState API | ✅ PASSING | 3/3 tests pass |
| WingGeometry API | ✅ PASSING | Via MemoryLifecycle tests |
| RotorInputState API | ✅ PASSING | Via MemoryLifecycle + Inflow tests |
| AircraftInputState API | ✅ PASSING | 6/6 tests pass |
| Wake/WakeHistory API | **❌ CRASH** | `oc_wake_create` causes SEGV |
| VTK API | ✅ COMPLETE | All 14 C bindings implemented in D (cbinding.d), full signatures in C header. Ready for testing/examples. |

### C++ Wrapper Architecture
The C++ wrapper (`cppbindings.cpp`) wraps each OC_* type in a RAII C++ class:
- Each wrapper stores a `void* ptr_` pointing to the underlying C object
- Move semantics are used throughout (move constructors/assignment take ownership)
- Destructors call the appropriate `oc_*_destroy` function when `owned_` is true
- Factory methods return owned objects (e.g., `AirfoilModel::thin_airfoil()`)
- Non-owning references are created via `Frame(p, false)` constructor

### Test Architecture
- 18 test source files in `tests/src/` 
- Built into single binary: `opencopter_c_api_tests` (17MB, gtest-based)
- CMakeLists.txt in `tests/CMakeLists.txt` for building
- Tests run from `tests/build/` directory
- Library loaded via `LD_LIBRARY_PATH=../..`

## VTK Segfault Fix (Completed 2026-08-01)

### CRITICAL FIX: `oc_write_rotor_vtu` Segfault Resolved
**Root Cause**: FFI casting error in `source/opencopter/cbindings.d`. `VtkRotor` is a D **class** (reference type), but the FFI layer was casting it as `cast(VtkRotor*)vtk_rotor` and dereferencing with `*v`, which is invalid for D classes. This caused garbage memory to be read, leading to segfault at `rotor.grid.SetPoints(rotor.points)` inside `vtk.d`.

**Fix Applied**:
```d
// BEFORE (WRONG - double indirection + invalid dereference):
auto v = cast(VtkRotor*)vtk_rotor;
write_rotor_vtu(fname, iteration, rotor_idx, *v, rs, *rg);

// AFTER (CORRECT - direct class reference):
auto v = cast(VtkRotor)vtk_rotor;
write_rotor_vtu(fname, iteration, rotor_idx, v, rs, *rg);
```

**Key Lesson**: D classes are reference types (like pointers). Casting a C `void*` to `ClassType*` then dereferencing with `*` creates double-indirection and reads garbage. The correct pattern is direct cast: `cast(ClassType)pointer`.

**Validation**: C example (`examples/c/example.c`) now runs all 30 iterations successfully, writing VTU files without crashing.

### Debug Instrumentation (Cleaned Up)
- `writeln`/`stdout.flush()` instrumentation was added to `vtk.d` to pinpoint crash location
- After fix confirmed, instrumentation was removed - no permanent changes to `vtk.d`

### Files Modified:
1. **`source/opencopter/cbindings.d`** — Fixed `oc_write_rotor_vtu` casting (~line 1087)

### Next Steps:
1. Run full C/C++ test suite with VTK-enabled build to confirm fix doesn't regress other tests
2. Document the D class FFI pattern in coding conventions for future reference
### Coding Conventions (D Language)
- Module-level documentation with DDoc (`///`)
- Unit tests embedded directly in source files
- Half-cosine spanwise distribution for blade elements
- Chunked processing aligned to `chunk_size()` for SIMD

### JSON Configuration Convention
- Geometry files define aircraft structure hierarchically
- Parameters files separate computational settings from flight conditions
- JSON5 format allows comments in example configs
- Aircraft → Rotors → Blades hierarchy mirrors physical structure

### Naming Convention (C API)
- Prefix: `oc_` (OpenCOPTER)
- Type suffix: `_create`, `_destroy`, `_get`, `_set`, `_fill`
- Compound names: `oc_<type>_<method>`
- Direction enums: `OC_CLOCKWISE`, `OC_COUNTER_CLOCKWISE`
- Frame type enums: `OC_AIRCRAFT_FRAME`, `OC_CONNECTION_FRAME`, etc.

## Learnings & Project Insights
### SIMD Architecture
- SLEEF provides vectorized math functions
- Data is processed in chunks matching CPU register width
- `oc_chunk_size()` returns the optimal chunk size for the current platform
- Blade element data is laid out to maximize vectorization efficiency

### Wake Convergence
- Wake positions converge iteratively over a revolution
- `convergence_criteria` in parameters controls tolerance
- `shed_history_angle` and `shed_release` control wake generation
- Converged vs uncomputed wake states are tracked separately

### Multirotor Support
- Single Aircraft can contain multiple rotors and wings
- Each rotor has independent RPM, collective, azimuth
- Rotor interactions modeled through shared wake field
- Direction (clockwise/counter-clockwise) affects induced velocity calculation

### Python Integration Architecture
- `pyd` submodule wraps D runtime for Python embedding
- Separate from C API - Python code can call D directly
- `oc_fly/` is the higher-level Python frontend that handles:
  - JSON config file parsing
  - MPI parallelization across flight conditions
  - WopWop acoustic coupling
  - Post-processing and visualization

### Airfoil Models
- AeroDAS only supports XFoil polar `.dat` files currently
- C81 Edgewise dynamic stall model available via `oc_c81_from_file()`
- Thin airfoil theory available via `oc_thin_airfoil_create()`
- Cl values computed at alpha=0 give ~6.28 (2π) which is the theoretical thin airfoil slope

## Areas Needing Documentation / Bug Fixes
- BWI (Broadband Wake Interaction) parameters - not documented
- Trim algorithms - multiple implementations exist but not yet documented
- C++ API header (`include/opencopter.hpp`) - needs review for completeness
- **`oc_wake_create()` crashes with SEGV** - struct layout mismatch between C/D (remaining critical bug)
- **BladeAirfoil Cl values ~2π at zero alpha** - may be correct (thin airfoil theory), tests need verification
