#! /bin/env bash

CPU_SUPPORT=native
BUILD_MODE=release

# Parse first argument: CPU target
if [ $# -gt 0 ]; then
    if [[ "$1" == "generic" || "$1" == "generic-avx" || "$1" == "generic-avx2" || "$1" == "generic-avx512f" ]]; then
        CPU_SUPPORT=$1
        if [ "$1" == "generic" ]; then
            CPU_SUPPORT=$1-avx
        fi
    elif [ "$1" == "native" ]; then
        CPU_SUPPORT=native
    elif [ "$1" == "--help" ]; then
        echo "Usage: ./build_linux.sh [cpu_target] [build_mode]"
        echo ""
        echo "CPU Target (default: native):"
        echo "  native            - Build for native CPU"
        echo "  generic           - Build for generic CPU (AVX)"
        echo "  generic-avx       - Build for generic CPU with AVX"
        echo "  generic-avx2      - Build for generic CPU with AVX2"
        echo "  generic-avx512f   - Build for generic CPU with AVX512F"
        echo ""
        echo "Build Mode (default: release):"
        echo "  debug             - Debug build (native CPU only)"
        echo "  release           - Release build"
        echo ""
        echo "Examples:"
        echo "  ./build_linux.sh                     # release + native"
        echo "  ./build_linux.sh native debug         # debug + native"
        echo "  ./build_linux.sh generic-avx2 release  # release + generic-avx2"
        exit 0
    else
        echo "Unrecognized CPU target argument: $1"
        echo "Use --help for usage information."
        exit 1
    fi
fi

# Parse second argument: build mode
if [ $# -gt 1 ]; then
    if [ "$2" == "debug" ]; then
        BUILD_MODE=debug
    elif [ "$2" == "release" ]; then
        BUILD_MODE=release
    else
        echo "Unrecognized build mode argument: $2"
        echo "Valid options: debug, release"
        echo "Use --help for usage information."
        exit 1
    fi

    # Check compatibility: debug builds only support native CPU
    if [ "$BUILD_MODE" == "debug" ] && [ "$CPU_SUPPORT" != "native" ]; then
        echo "Error: Debug builds only support native CPU targets."
        echo "Please use 'native' with debug mode, or switch to release mode for generic CPU builds."
        exit 1
    fi
fi

# Check for anaconda
if command -v conda &> /dev/null; then
    echo "Anaconda/Miniconda is installed."
else
    echo "Anaconda/Miniconda is not installed. Please install Anaconda/Miniconda before continuing"
    exit -1
fi

# Check for opencopter environment existence
if conda list | grep opencopter; then
    echo "opencopter environment already exists, skipping creation."
else
    conda env create -f environment.yml
fi

# Check for opencopter environment activation
if [ "$CONDA_DEFAULT_ENV" == "opencopter" ]; then
    echo "opencopter environment already activated, skipping activation."
else
    conda activate opencopter
fi

# Configure build type
BUILD_CONFIG=$BUILD_MODE-$CPU_SUPPORT
if [ "$CPU_SUPPORT" == "native" ]; then
    if grep avx512 /proc/cpuinfo &> /dev/null; then
        echo "AVX512 support enabled."
        BUILD_CONFIG=$BUILD_CONFIG-512
    fi
fi

echo "Building OpenCOPTER with build configuration $BUILD_CONFIG"

dub build -c library-python312 -b $BUILD_CONFIG --compiler=ldc2

cd oc_fly/dependencies/wopwopd

dub build -c library-python312 -b $BUILD_CONFIG --compiler=ldc2

cd -
