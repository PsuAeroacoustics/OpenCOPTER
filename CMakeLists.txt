cmake_minimum_required(VERSION 3.10)
project(OpenCOPTER CXX)

# BUILD_TYPE comes from dub and carries dub's build type verbatim, e.g.
# "release-native", hence the substring matches rather than equality.
if("${BUILD_TYPE}" STREQUAL "")
	set(CMAKE_BUILD_TYPE "Debug")
elseif("${BUILD_TYPE}" MATCHES "debug")
	set(CMAKE_BUILD_TYPE "Debug")
elseif("${BUILD_TYPE}" MATCHES "release")
	set(CMAKE_BUILD_TYPE "Release")
endif()

add_library(cppbindings OBJECT ${CMAKE_SOURCE_DIR}/source/opencopter/cppbindings.cpp)
target_compile_options(cppbindings PRIVATE -fPIE -fPIC)
target_compile_features(cppbindings PUBLIC cxx_std_23)
target_include_directories(cppbindings PRIVATE ${CMAKE_SOURCE_DIR}/include)
