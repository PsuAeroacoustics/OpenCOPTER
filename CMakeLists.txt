cmake_minimum_required(VERSION 3.10)
project(OpenCOPTER CXX)

if(${BUILD_TYPE} STREQUAL "")
	set(CMAKE_BUILD_TYPE "Debug")
elseif(${BUILD_TYPE} MATCHES "debug")
	set(CMAKE_BUILD_TYPE "Debug")
elseif(${BUILD_TYPE} MATCHES "release")
	set(CMAKE_BUILD_TYPE "Release")
endif()

set(CMAKE_CXX_STANDARD 23)

add_custom_target(cppbindings ALL
    COMMAND ${CMAKE_CXX_COMPILER} -c -fPIC --std=c++23
            -I${CMAKE_SOURCE_DIR}/include
            ${CMAKE_SOURCE_DIR}/source/opencopter/cppbindings.cpp
            -o ${CMAKE_SOURCE_DIR}/source/opencopter/cppbindings.o
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMENT "Compiling cppbindings.cpp")