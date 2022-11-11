cmake_minimum_required(VERSION 3.9)

project(Terrain_Trees_core)
set(terrain_trees_VERSION_MAJOR 1)
set(terrain_trees_VERSION_MINOR 0)

SET(CMAKE_MODULE_PATH ${PROJECT_BINARY_DIR}/cmake/Modules/)
MESSAGE(STATUS "local cmake root: " ${CMAKE_MODULE_PATH})

set(LANGUAGE C++)
#set(CMAKE_VERBOSE_MAKEFILE on)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib )

set(CMAKE_BUILD_TYPE Release)

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
#find_package(Doxygen REQUIRED)
find_package(OpenMP)




#if(CMAKE_COMPILER_IS_GNUCXX) #do not work on mac
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "-O3 -std=c++11 -fopenmp") ## Optimize
#    set(CMAKE_EXE_LINKER_FLAGS "-s -lrt") ## Strip binary
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(CMAKE_CXX_FLAGS "-O3 -march=native -std=c++11 -fopenmp") ## Optimize
    set(CMAKE_EXE_LINKER_FLAGS "-s") ## Strip binary
endif()

include_directories(sources)
add_subdirectory(sources)

target_include_directories(_terrain_trees_core PUBLIC ${PROJECT_BINARY_DIR}/sources/)
target_include_directories(_terrain_trees_apps PUBLIC ${PROJECT_BINARY_DIR}/sources/)
