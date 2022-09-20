# OpenCOPTER
OpenCOPTER (COupled Potential Theory Extensions for Rotors) is a library for fast and efficient simulation of multirotor aerodynamics.

## Dependencies


OpenCOPTER has a number of dependencies that need to be setup and installed first. Some of these are provided as git submodules and are largely transparent to the user. If you are cloning a fresh repository run:
```
	git clone --recurse-submodules https://github.com/PsuAeroacoustics/OpenCOPTER
```

If the repository is already cloned, but submodules were not cloned, run the following command to acquire the submodules:

```
	git submodule update --init
```

Setting up the remaining dependencies is described below.

### Required Dependencies

- The LLVM D Compiler (ldc) (https://github.com/ldc-developers/ldc)
- SLEEF (https://sleef.org/)
- OpenBLAS (https://www.openblas.net/)
- Lapack
- Lapacke

**Linux (Ubuntu)**

Most of these dependencies can be obtained through the apt package manager on Ubuntu. Make sure the development version of the packages are installed:

```
	sudo apt install ldc
	sudo apt install libopenblas-dev
	sudo apt install liblapack-dev
	sudo apt install liblapacke-dev
```

The exception to this is SLEEF. As OpenCOPTER statically links SLEEF, when configuring it with cmake, make sure that the `-DBUILD_SHARED_LIBS=OFF` flag is passed to cmake. OpenCOPTER expects that SLEEF is installed in `/usr/local/lib` which is the default prefix on linux and macOS

**macOS**

Homebrew can be used to install most of these dependencies. The exception is SLEEF, which while there is a brew formula, it doesn't build the static library. This can be built the same way as the linux version. It is expected to be installed in `/usr/local/lib` on macOS as well.

Homebrew can also be used to install openblas. On macOS openblas includes lapack and lapacke, so no further dependencies need to be installed. If Homebrew is not an option, openblas can be built from source and installed to `/usr/local/lib`.

Using Homebrew:
```
	brew install dub
	brew install ldc
	brew install -s openblas
```

The `-s` option for the openblas install makes Homebrew build the library from source, which is required to avoid openmp runtime errors.

**Windows**

Unfortunately OpenCOPTER does not support Windows natively at this time. Work is being done to get support going but for the time being, WSL and WSL2 can be used to run this code using the linux directions above.

### Optional Dependencies

- cmake 3.14 or greater
- libvtk-9.1

If building the library with vtk support, we currently explicitly link with libvtk-9.1 as this is what is offered in the Ubuntu 22.04 repositories. The libraries are unfortunately installed with an explicit `-9.1` suffix and no symlinks, which leaves us linking to a specific version. At the moment, the easiest way to install this version is through the git repository as the pre-built packages on the kitware site are for 9.2 at the time of this writing.

### Example/Validation Project Dependencies


The example project and boxwell validation project have a few more dependencies that need to be put in place if the user is interested in running them. The example `example/python/example.py` expects the following python packages to be installed:

- matplotlib
- numpy

The boxwell validation project expects two more dependencies to be installed. The first is our acoustics wopwop coupling library, wopwopd (https://github.com/PsuAeroacoustics/wopwopd.git). This can be placed anywhere but must be registered to the local dub repository. This can be done by navigating to the repository on the command line and running the following command:

```
	dub add-local ./
```
The second dependency is our fork of the D language binding for matplotlib, matplotlib-d (https://github.com/Rob-Rau/matplotlib-d). The validation project expects this to be installed in the same directory as the OpenCOPTER directory. This package requires that the python package matplotlib be installed.

No other action needs to be taken with these dependencies as their build process is taken care of when building the validation project.

## Building

Once all the dependencies have been installed, run the following command to build the library:

```
	dub build -c library-novtk -b release --compiler=ldc2
```

This command builds the library without vtk support as a portable release build using the ldc compiler. There are a number of different configurations and build types that can be used.

### Build Types (`-b`)


| Configuration                      | Description                                    |
|------------------------------------|------------------------------------------------|
| debug                              | Basic debug mode, no optimizations.            |
| release                            | Portable release build. Non-CPU specific optimizations enabled |
| release-native                     | Native release build. Optimizes for host CPU   |
| release-native-512                 | Native release build. Optimizes for host CPU, Encourages use of AVX512 instructions, suppressing compiler instruction heuristics |
| release-generic-avx                | Portable release mode with AVX enabled         |
| release-generic-avx2               | Portable release mode with AVX2 enabled        |
| release-generic-avx512f            | Portable release mode with AVX512F enabled     |

OpenCOPTER has been designed so that the compilers auto-vectorizer can be judiciously employed. This means that there may be large performance gains by using the `-native` build types. However it is important to note that these builds will *not* be portable.

### Configurations (`-c`)

| Build Type                    | Description                                                                 |
|-------------------------------|-----------------------------------------------------------------------------|
| library                       | Builds the dynamic library for use with other D code. No python wrappers included |
| library-novtk                 | Builds the dynamic library for use with other D code. Python wrappers and VTK support excluded |
| library-python<version>       | Builds the dynamic library for use with other D code. Python wrappers are built for the specific version of python. <version> can be any of 33, 34, 35, 36, 37, 38, 39, or 310 |
| library-python<version>-novtk | Builds the dynamic library for use with other D code. VTK support excluded. Python wrappers are built for the specific version of python. <version> can be any of 33, 34, 35, 36, 37, 38, 39, or 310 |

## Running the Examples


All that needs to be done to run the example after libopencopter has been built with python support is running the following command in the examples/python directory:

```
	python3 example.py
```

The HART-II validation can also be referenced as an example project. It has a number of command line arguments that can be set to control the collective pitch of the rotor, the lateral and longitudinal cyclic, and more. To get a full list of command line options, on the command line run:

```
	./hart_val -h
```
