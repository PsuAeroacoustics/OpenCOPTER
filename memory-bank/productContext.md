# Product Context: OpenCOPTER

## Why This Project Exists
Aircraft rotor aerodynamics simulation traditionally requires computationally expensive CFD for accurate results. OpenCOPTER addresses this by providing a fast vortex-lattice/blade-element method (VLM/BEM) that achieves good accuracy at a fraction of the computational cost, enabling rapid design space exploration for rotary-wing and eVTOL aircraft configurations.

## Problems Solved
1. **Design Cycle Speed**: Traditional CFD takes hours/days per configuration. OpenCOPTER enables minute-level iterations during concept design.
2. **Multirotor Interaction**: Accurately models rotor-rotor and rotor-wing interference effects critical to eVTOL/tiltrotor designs.
3. **Cross-Language Access**: Provides unified C/C++/Python APIs so engineers can use their preferred language without losing functionality.
4. **Aeroacoustic Prediction**: Couples aerodynamics with wopwop for acoustic predictions (discrete frequency and broadband noise).

## How It Works (User Perspective)
1. **Define geometry**: JSON file specifies rotors, blades, wings with radii, origins, chord distributions, twist, airfoils.
2. **Define flight conditions**: JSON/JSON5 file specifies collective, RPM, airspeed, density, angle of attack.
3. **Run simulation**: `oc_fly` processes input files and calls OpenCOPTER's C library for computation.
4. **Post-process**: Results include thrust/torque coefficients, blade loads, wake visualization (.vtu), and optionally acoustic predictions via wopwopd.

## User Experience Goals
- **Minimal setup**: Anaconda environment + build script should get users started in minutes.
- **Intuitive input**: JSON-based configuration files mirror physical aircraft terminology (rotors, blades, wings).
- **Clear validation**: HART-II and other benchmark cases included for verification.
- **Visualization**: VTK output integrable with Paraview for wake/aerodynamic inspection.

## Target Aircraft Types
- Single/dual rotors (helicopter mode)
- Tiltrotor aircraft (V22 Osprey class)
- Tiltwing aircraft
- Multi-rotor eVTOL (vertical takeoff and landing)
- Lifting rotor + wing configurations
- Notional concept designs

## Key Performance Characteristics
- **Parallelizable**: oc_fly supports MPI-based parallel execution across flight conditions.
- **AVX-optimized**: SIMD-accelerated math routines via SLEEF library.
- **Chunked processing**: Vectorized operations aligned with CPU cache lines.