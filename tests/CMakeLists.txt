cmake_minimum_required(VERSION 3.14)
project(opencopter_c_api_tests C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Options
# OFF by default: libopencopter is built by dub without ASan, so enabling this
# without also instrumenting the library breaks linking. See the note below.
option(ENABLE_ASAN "Enable AddressSanitizer for memory error detection" OFF)
option(ENABLE_VALGRIND_TESTS "Enable Valgrind memcheck tests" ON)

# Find GoogleTest via FetchContent (downloads if not already present)
include(FetchContent)
FetchContent_Declare(
  googleTest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG v1.14.0
  GIT_SHALLOW TRUE
)
# Prevent GoogleTest from overriding our compiler/linker options
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googleTest)

# Locate the OpenCOPTER shared library and headers
# The library is produced by `dub build` and lives in the project root
set(OPENCOPTER_ROOT "${CMAKE_SOURCE_DIR}/..")

find_library(OPENCOPTER_LIB NAMES opencopter libopencopter.so
             PATHS "${OPENCOPTER_ROOT}"
             NO_DEFAULT_PATH)

if(NOT OPENCOPTER_LIB)
  message(WARNING "OpenCOPTER library not found at ${OPENCOPTER_ROOT}. "
                  "Make sure you have built OpenCOPTER first (e.g. ./build_linux.sh). "
                  "You can override the search path by setting OPENCOPTER_LIB manually.")
endif()

include_directories("${OPENCOPTER_ROOT}/include")

# Enable testing
enable_testing()

# -----------------------------------------------------------------------
# Test source files – each focuses on memory safety of a slice of the API
# -----------------------------------------------------------------------
set(TEST_SOURCES
  # Original memory-safety tests
  src/test_memory_lifecycle.cpp
  src/test_memory_frame.cpp
  src/test_memory_aircraft.cpp
  src/test_memory_blade.cpp
  src/test_memory_wake.cpp
  src/test_memory_buffer.cpp
  # Phase 2: Config, Airfoil, BladeAirfoil
  src/test_api_config.cpp
  src/test_api_airfoil.cpp
  src/test_api_bladeairfoil.cpp
  # Phase 3: BladeGeometry chain
  src/test_api_bladegroup.cpp
  src/test_api_rotorgeo_advanced.cpp
  # Phase 4: Inflow & State chain
  src/test_api_inflow.cpp
  src/test_api_aircraftstate.cpp
  # Phase 5: Wake & Simulation
  src/test_api_wake.cpp
  src/test_api_simulation.cpp
  # Phase 6: Utility & VTK
  src/test_api_utility.cpp
  src/test_api_vtk.cpp
  # Phase 7: Extended sentinel tests
  src/test_api_sentinel.cpp
  # Phase 8: Children API tests
  src/test_api_children.cpp
  # Phase 9: C++ RAII wrapper tests
  src/test_cpp_api.cpp
)

# Create single test executable
add_executable(opencopter_c_api_tests ${TEST_SOURCES})

target_link_libraries(opencopter_c_api_tests
  PRIVATE
    GTest::gtest
    GTest::gtest_main
    pthread
    m
)

# Link against the OpenCOPTER shared library when available
if(OPENCOPTER_LIB)
  target_link_libraries(opencopter_c_api_tests PRIVATE ${OPENCOPTER_LIB})
  set_target_properties(opencopter_c_api_tests PROPERTIES
    BUILD_RPATH "${OPENCOPTER_ROOT}"
    INSTALL_RPATH "${OPENCOPTER_ROOT}")
endif()

# AddressSanitizer flags (compile + link).
# NOTE: The OpenCOPTER shared library is built by `dub` without ASan, so
# enabling ASan here will cause conflicts.  Keep it OFF by default and
# only turn it on when the host library is also rebuilt with -fsanitize=address.
if(ENABLE_ASAN)
  message(WARNING "ASan enabled — only meaningful if libopencopter.so was also compiled with -fsanitize=address")
  target_compile_options(opencopter_c_api_tests PRIVATE -fsanitize=address -fno-omit-frame-pointer)
  target_link_options(opencopter_c_api_tests PRIVATE -fsanitize=address)
endif()

# Register with CTest so `ctest` discovers individual GoogleTest cases
include(GoogleTest)
gtest_discover_tests(opencopter_c_api_tests DISCOVERY_TIMEOUT 120)

# Optional: Valgrind memcheck test
if(ENABLE_VALGRIND_TESTS)
  find_program(MEMCHECK_VALGRIND valgrind)
  if(MEMCHECK_VALGRIND)
    add_test(
      NAME memcheck_full_suite
      COMMAND ${MEMCHECK_VALGRIND}
              --tool=memcheck
              --leak-check=full
              --errors-for-leak-kinds=definite,possible
              --error-exitcode=1
              $<TARGET_FILE:opencopter_c_api_tests>
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    )
  endif()
endif()