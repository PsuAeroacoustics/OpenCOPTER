# System Patterns: OpenCOPTER

## Architecture Overview
OpenCOPTER follows a layered architecture:
```
Input (JSON/JSON5) → oc_fly (Python frontend) → C API → D Core Library
                                                                        ├─ Wake Dynamics
                                                                        ├─ Blade Element Analysis
                                                                        ├─ Inflow Models
                                                                        └─ Aerodynamic Forces
                                                                    ↓
                                                              VTK Output / Acoustics
```

## Core Components

### 1. Aircraft Model (`aircraft/`)
- **geometry.d**: Hierarchical frame tree (root → rotors/wings → blades)
- **input.d**: Flight condition inputs (collective, RPM, cyclic, freestream)
- **state.d**: Runtime state (C_T, C_Q, velocities, angles of attack per blade element)

### 2. Wake System (`wake.d`, `bladeelement.d`)
- Lagrangian vortex filament tracking
- Shed and release wake criteria configurable per rotor
- Half-cosine spanwise distribution for blade elements
- Chunked memory layout for SIMD vectorization

### 3. Inflow Models (`inflow/`)
- **huangpeters/**: Full Beddoes-Leishman dynamic stall model with system/iteration solvers
- **beddos.d**: Simplified Beddoes-Leishman model
- **wingingflow.d**: Wing-induced velocity calculations
- **simplewing.d**: Simplified wing inflow for fast simulation

### 4. Airfoil Models (`airfoilmodels/`)
- **aerodas.d**: AeroDAS lookup table interpolation (Cl, Cd vs alpha, Mach)
- **c81.d**: Edgewise dynamic stall model (C81 polar format)
- **thinaf.d**: Thin airfoil theory (constant Cl_alpha, zero drag)

### 5. Lift Models (`liftmodels/`)
- **steadylift.d**: Steady lift surface calculations
- **vortexlattice.d**: Vortex lattice method for wing aerodynamics
- **weissingerl.d**: Weissinger L-method for low-speed lifting surface theory

### 6. Math Utilities (`math/`)
- **sleef.d**: SLEEF SIMD math functions (trig, sqrt, etc.)
- **blas.d**: BLAS routines via OpenBLAS
- **lapacke.d**: LAPACK linear algebra solver bindings
- **vectorarray.d**: Vector array operations for blade elements

### 7. VTK Output (`vtk.d`)
- .vtu file generation for rotors, wings, wakes
- Wake field visualization (induced velocity, vorticity)
- Per-frame and per-simulation output modes

### 8. Trim (`trim.d`)
- Collective trim to match target C_T per rotor
- Iterative correction with convergence criteria

## Design Patterns Used

### Opaque Pointer / PIMPL Pattern
All complex types (OC_Aircraft, OC_Wake, OC_Inflow) use opaque pointers in the C API. Implementation details are hidden behind the ABI boundary.

### Factory Functions
Objects are created via factory functions (e.g., `oc_huang_peters_create`, `oc_aero_das_from_xfoil_polar`) rather than direct struct construction.

### Strategy Pattern for Inflow Models
`OC_Inflow` is an abstract interface implemented by different strategies:
- `OC_HuangPeters` - full dynamic stall
- `OC_NullInflow` - zero induced velocity
- `OC_WingInflow` - wing-induced flow

### Chunked SIMD Processing
Arrays are processed in chunks (`oc_chunk_size()`) aligned with CPU vector register widths. SLEEF provides vectorized math functions operating on these chunks.

### Frame Hierarchy
Frames form a tree structure for coordinate transforms:
```
Aircraft Root Frame
├── Rotor 1 Frame → Rotor 1 Geometry → Blades [0..N]
├── Rotor 2 Frame → Rotor 2 Geometry → Blades [0..N]
└── Wing 1 Frame → Wing Geometry → Parts [0..N]
```

## Data Flow

### Simulation Step (`oc_simulation_step`)
1. Update inflow models with current wake state
2. Compute induced velocities at blade element points
3. Calculate angle of attack and effective velocity for each element
4. Look up airfoil coefficients (Cl, Cd) from polar data
5. Integrate blade element forces (dC_T, dC_Q, dC_L, dC_D)
6. Update wake filaments based on shed/release criteria
7. Apply convergence check on wake positions

### Wake Update Flow
```
oc_simulation_step → oc_inflow_update
    → compute_velocities_at_blade_points
    → update_vortex_filament_positions
    → converge_wake (iteration loop)
```

## Key Technical Decisions

### D Language Choice
- Garbage collection for complex data structures
- Built-in unit testing (`unittest` blocks)
- Compile-time metaprogramming for code generation
- Foreign function interface to C is seamless

### SIMD via SLEEF over Intrinsics
SLEEF provides auto-vectorizing functions that work across compilers without platform-specific intrinsics. This enables portable vectorization while maintaining performance.

### JSON over XML/JSON5 for Config
oc_fly uses JSON5 (with comments) for human-editable config files. Programmatic interfaces use standard JSON.

### VTK as Primary Visualization
VTK .vtu format provides:
- Cross-platform visualization (ParaView, VisIt)
- Native Python library support (PyVista)
- Standard CFD/CAE interchange format