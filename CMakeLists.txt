# NodeGenLib - A library for advancing front node generation
# Copyright (C) 2024 Bruno Degli Esposti
#
# This file is part of NodeGenLib.
#
# NodeGenLib is licensed under the GNU Lesser General Public License (LGPL), 
# either version 3 or (at your option) any later version. See the LICENSE 
# file in the root of the project for more details.
#
# This file is distributed without any warranty.

# https://cmake.org/cmake/help/latest/manual/cmake-commands.7.html
cmake_minimum_required(VERSION 3.24 FATAL_ERROR)

project("NodeGenLib" CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
add_compile_options(-Wall -Wextra -Wpedantic)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

include(FetchContent)
FetchContent_Declare(
  nanoflann
  GIT_REPOSITORY https://github.com/jlblancoc/nanoflann.git
  GIT_TAG origin/master
  FIND_PACKAGE_ARGS NAMES nanoflann
)
FetchContent_MakeAvailable(nanoflann)

find_package(Matlab REQUIRED)
find_package(Eigen3 3.4 REQUIRED)
find_package(OpenCASCADE REQUIRED)
include_directories(SYSTEM ${OpenCASCADE_INCLUDE_DIR})
link_directories(${OpenCASCADE_LIBRARY_DIR})
if(UNIX AND NOT APPLE)
    set(CMAKE_BUILD_RPATH ${OpenCASCADE_LIBRARY_DIR})
endif()

matlab_add_mex(
	NAME advancing_front
	SRC src/advancing_front.cpp
	R2018a
)

matlab_add_mex(
    NAME step_generate_boundary_nodes
    SRC src/step_generate_boundary_nodes.cpp
    R2018a
)

matlab_add_mex(
    NAME step_plot_pcurves
    SRC src/step_plot_pcurves.cpp
    R2018a
)

matlab_add_mex(
    NAME step_print_contents
    SRC src/step_print_contents.cpp
    R2018a
)

matlab_add_mex(
    NAME step_print_faces_info
    SRC src/step_print_faces_info.cpp
    R2018a
)

target_link_libraries(advancing_front Eigen3::Eigen)
target_link_libraries(advancing_front nanoflann::nanoflann)
target_link_libraries(step_generate_boundary_nodes Eigen3::Eigen)
target_link_libraries(step_generate_boundary_nodes nanoflann::nanoflann)
foreach(LIB ${OpenCASCADE_LIBRARIES})
	target_link_libraries(step_generate_boundary_nodes ${LIB})
	target_link_libraries(step_plot_pcurves ${LIB})
	target_link_libraries(step_print_contents ${LIB})
    target_link_libraries(step_print_faces_info ${LIB})
endforeach()
