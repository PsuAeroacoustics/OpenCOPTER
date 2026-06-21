#! /bin/env bash

CPU_SUPPORT=native
if [ $# -gt 0 ]; then
    if [[ "$1" == "generic" || "$1" == "generic-avx" || "$1" == "generic-avx2" || "$1" == "generic-avx512f" ]]; then
        CPU_SUPPORT=$1
        if [ "$1" == "generic" ]; then
            CPU_SUPPORT=$1-avx
        fi
    elif [ "$1" == "native" ]; then
        CPU_SUPPORT=native
    elif [ "$1" == "--help" ]; then
        echo "Usage:"
        echo "  For native CPU build:"
        echo "      ./build_linux.sh"
        echo "          or"
        echo "      ./build_linux.sh native"
        echo "  For generic CPU build:"
        echo "      ./build_linux.sh generic"
        echo "  For generic CPU build with AVX support:"
        echo "      ./build_linux.sh generic-avx"
        echo "  For generic CPU build with AVX2 support:"
        echo "      ./build_linux.sh generic-avx2"
        echo "  For generic CPU build with AVX512F support:"
        echo "      ./build_linux.sh generic-avx512f"
        exit 0
    else
        echo "Unrecognized input argument: $1"
        exit -1
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
BUILD_CONFIG=release-$CPU_SUPPORT
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
