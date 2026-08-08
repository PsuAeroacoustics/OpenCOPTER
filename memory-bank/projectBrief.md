# Project Brief: OpenCOPTER

## Overview
**OpenCOPTER** (COupled Potential Theory Extensions for Rotors) is a high-performance library for fast and efficient multirotor aerodynamics simulation. It uses vortex lattice methods and blade element theory to model rotor/wing interactions with high computational efficiency.

## Core Requirements
- Simulate aerodynamics of multi-rotor aircraft configurations (tiltrotor, tiltwing, eVTOL, etc.)
- Provide coupled wake/geometry aeroacoustic predictions
- Support multiple airfoil models including dynamic stall (Edgewise Stall model)
- Generate VTK visualization output (.vtu files) for computational fluid dynamics post-processing
- Expose C, C++, and Python APIs for cross-language accessibility

## Key Capabilities
1. **High-Resolution Wake Dynamics**: Lagrangian vortex wake with fine timestep formulation specifically designed for acoustics-grade airload resolution
2. **Blade Element Theory**: Spanwise discretized blade aerodynamics with airfoil polars (XFoil format)
3. **Inflow Models**: Huang-Peters (Beddoes-Leishman), and wing inflow models
4. **Wing Modeling**: Weissinger L lifting surface method for wing aerodynamics
5. **VTK Output**: VTK XML Unstructured Grid (.vtu) output for rotors, wings, wakes, and wake fields
6. **Trim Capability**: Collective/cyclic trim to match target thrust coefficients
7. **BWIAcoustics**: Broadband and discrete noise predictions

## Validation Targets
- HART II rotor cases (hover and forward flight)
- Hart II single rotor configurations
- Experimental validation against published datasets

## Project Structure
```
OpenCOPTER/
├── source/opencopter/        # Core D library source
│   ├── aircraft/              # Aircraft geometry, input, state
│   ├── airfoilmodels/         # Airfoil models (AeroDAS, C81, Thin Airfoil)
│   ├── inflow/                # Inflow models (Huang-Peters, WingInflow, etc.)
│   ├── liftmodels/            # Lift surface methods (Weissinger L)
│   └── math/                  # Math utilities (BLAS, LAPACK, Sleef SIMD)
├── include/                   # C/C++ public headers
├── examples/                  # Usage examples (C, C++, Python)
├── tests/                     # Test suite
├── oc_fly/                    # Python frontend for batch simulation
├── dependencies/              # Git submodules (numd, pyd, vtkd)
└── validation/                # Validation cases and reference data
```

## License
MIT License

## Repository
- **Remote**: `git@github.com:PsuAeroacoustics/OpenCOPTER.git`
- **Documentation**: https://psuaeroacoustics.github.io/