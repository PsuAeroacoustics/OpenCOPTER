cmake_minimum_required(VERSION 3.10)
project(OpenCOPTER CXX)

set(CMAKE_CXX_STANDARD 23)

add_custom_target(cppbindings ALL
    COMMAND ${CMAKE_CXX_COMPILER} -c -fPIC
            -I${CMAKE_SOURCE_DIR}/include
            ${CMAKE_SOURCE_DIR}/source/opencopter/cppbindings.cpp
            -o ${CMAKE_SOURCE_DIR}/source/opencopter/cppbindings.o
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMENT "Compiling cppbindings.cpp")