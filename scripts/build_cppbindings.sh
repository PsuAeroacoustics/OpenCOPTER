#! /bin/sh

# Builds the C++ bindings object that dub links into libopencopter.
#
# Run from dub's preBuildCommands. This has to fail loudly: a silent failure
# here leaves a stale cppbindings.cpp.o in build/ that dub will happily link,
# which looks like a successful build of the wrong sources.

set -e

cd "$(dirname "$0")/.."

BUILD_DIR=build

# Prefer Ninja when it is installed, otherwise fall back to CMake's default
# generator, which is always available wherever CMake is.
if command -v ninja > /dev/null 2>&1; then
    GENERATOR=Ninja
else
    GENERATOR="Unix Makefiles"
fi

# dub exports its build type as DUB_BUILD_TYPE. The plain BUILD_TYPE this
# script used to read was never set by dub, so CMakeLists.txt always saw an
# empty value and its if() failed to parse.
BUILD_TYPE=${DUB_BUILD_TYPE:-${BUILD_TYPE:-release}}

configure() {
    cmake . -DBUILD_TYPE="$BUILD_TYPE" -B "$BUILD_DIR" -G "$GENERATOR" -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
}

mkdir -p "$BUILD_DIR"

# CMake refuses to configure over a cache written by a different generator.
# That is recoverable -- everything in build/ is generated -- so start over.
if ! configure; then
    echo "cmake configure failed, reconfiguring $BUILD_DIR from scratch" >&2
    rm -rf "$BUILD_DIR"
    configure
fi

cmake --build "$BUILD_DIR"
