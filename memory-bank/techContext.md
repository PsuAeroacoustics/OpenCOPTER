# Technical Context: OpenCOPTER

## Language & Compiler
- **Primary language**: D (DLang) - source code in `.d` files
- **Compiler**: LDC (LLVM-based D compiler) - `ldc2`
- **C++ binding file**: `cppbindings.cpp` (compiled via CMake)
- **C++ standard**: C++23 (`cxx_std_23`)
- **Build system**: Dub (D package manager) + CMake (for C++ bindings/VTK)

## Dependencies

### Core D Dependencies (git submodules in `dependencies/`)
| Dependency | Purpose |
|------------|---------|
| `numd` | Numerical computing library for D (BLAS/LAPACK wrappers, math utilities) |
| `pyd` | Python-D bindings for embedding OpenCOPTER in Python |
| `vtkd` | VTK bindings for D (visualization output) |
| `kxml` | XML parsing for D language |

### System Dependencies (Linux)
| Package | Purpose |
|---------|---------|
| `libopenblas-dev` | BLAS linear algebra routines |
| `liblapack-dev` | LAPACK linear algebra solver |
| `liblapacke-dev` | LAPACK C interface |
| `libsleef-dev` | SLEEF SIMD math library (trig, sqrt, etc.) |
| `libvtk-9.1` | VTK visualization (optional, `-vtk` configs) |
| `cmake ≥ 3.10` | CMake build system for cppbindings and vtkd |

### Python Dependencies (oc_fly)
| Package | Purpose |
|---------|---------|
| `numpy` | Array operations, math |
| `scipy` | Scientific computing |
| `mpi4py` | MPI parallelization |
| `pyjson5` | JSON5 parsing (comments in config files) |

### Validation Project Dependencies
| Dependency | Purpose |
|------------|---------|
| `wopwopd` | Acoustics coupling library (github.com/PsuAeroacoustics/wopwopd) |
| `matplotlib-d` | Plotting (Rob-Rau/matplotlib-d fork) |

## Build System

### Dub Configurations
| Configuration | Description |
|---------------|-------------|
| `library` | Dynamic library + VTK support, no Python wrappers |
| `library-novtk` | Dynamic library without VTK |
| `library-python3X` | Python bindings for specific Python version (33-314) |
| `library-python3X-novtk` | Python bindings without VTK |
| `library_huang` | Huang-Peters inflow variant |

### Build Types
| Type | Options | Use Case |
|------|---------|----------|
| `debug` | debugMode, debugInfo | Development debugging |
| `debug-native` | + AddressSanitizer | Memory error detection |
| `release` | optimize, inline | Portable production build |
| `release-native` | -mcpu=native | Host-optimized (not portable) |
| `release-generic-avx` | AVX enabled | Portable SIMD acceleration |
| `release-generic-avx2` | AVX2 enabled | Portable AVX2 acceleration |
| `release-generic-avx512f` | AVX512F enabled | AVX512 portable build |
| `release-native-512` | EVEX512, AMX, no boundscheck | Maximum performance (Linux only) |

### Build Command
```bash
# Standard build
dub build -c library-novtk -b release --compiler=ldc2

# With anaconda + build script
./build_linux.sh native
```

## Platform Support
| Platform | Status | Notes |
|----------|--------|-------|
| Linux (Ubuntu 22.04+) | ✅ Full | Primary development platform |
| macOS | ✅ Partial | Accelerate framework, Homebrew dependencies |
| Windows | ❌ Not native | Requires WSL/WSL2 |

## Python Integration
- Via `pyd` submodule - generates `.so` with D runtime embedded
- Python version-specific configurations (3.3 through 3.14)
- After build: `libopencopter.so` can be imported directly in D code that links to it
- `oc_fly/` is a separate Python package that wraps the C library via ctypes/ffi

## Memory Management
- D garbage collection for high-level types (Aircraft, Wake, etc.)
- Manual memory management through C API (`destroy` functions)
- Chunk-based arrays for blade elements aligned to SIMD width
- `memory.d` provides custom allocation utilities

## Key File Paths
```
source/
├── opencopter/
│   ├── cbindings.d          # C FFI bindings (entry point)
│   ├── cppbindings.cpp      # C++ bindings (compiled via CMake)
│   ├── wake.d               # Vortex wake system
│   ├── bladeelement.d       # Blade element analysis
│   ├── trim.d               # Thrift/collective trim
│   ├── vtk.d                # VTK output
│   ├── atmosphere.d         # Atmospheric properties
│   ├── bwi.d                # Broadband noise
│   └── aircraft/            # Aircraft model components
├── dependencies/            # Git submodules
include/                     # C public API header
tests/                       # D unit tests + C++ tests
oc_fly/                      # Python frontend
validation/                  # Validation cases
```

## Environment Variables
| Variable | Purpose |
|----------|---------|
| `BUILD_TYPE` | Passed to CMake for build type selection |
| `$CONDA_PREFIX` | Searched for libraries in anaconda environments |
| `$HOME/.local/lib` | Searched for user-installed libraries |

## Testing Framework
- **D tests**: Built-in `unittest` blocks in source files
- **C++ tests**: Separate test suite in `tests/` using CMake
- **Test frameworks used**:
  - D native `unittest` blocks
  - Google Test (gtest) for C++ API tests