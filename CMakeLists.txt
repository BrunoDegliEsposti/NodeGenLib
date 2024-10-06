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
cmake_minimum_required(VERSION 3.12)

project(NodeGenLib VERSION 0.2 LANGUAGES CXX)

find_package(Matlab REQUIRED)
find_package(Eigen3 3.4 REQUIRED)
find_package(OpenCASCADE REQUIRED)

option(NANOFLANN_BUILD_EXAMPLES "Build nanoflann examples" OFF)
option(NANOFLANN_BUILD_TESTS "Build nanoflann unit tests" OFF)
include(FetchContent)
FetchContent_Declare(
  nanoflann
  GIT_REPOSITORY https://github.com/jlblancoc/nanoflann.git
  GIT_TAG v1.6.1
)
FetchContent_MakeAvailable(nanoflann)

# MATLAB contains its own copy of OpenCASCADE and when our MEX
# functions are called on linux, MATLAB's version may get
# dynamically loaded instead of the system-wide version,
# because MATLAB puts its own directories at the start of LD_LIBRARY_PATH.
# The problem is solved by explicitly setting the RPATH in the MEX files.
if(UNIX AND NOT APPLE)
    set(CMAKE_BUILD_RPATH "${OpenCASCADE_LIBRARY_DIR}")
endif()

# Build instructions for MEX function advancing_front
matlab_add_mex(
	NAME advancing_front
	SRC src/advancing_front.cpp
	R2018a
)
target_compile_features(advancing_front
    PUBLIC cxx_std_17)
target_compile_options(advancing_front
    PUBLIC -Wall -Wextra -Wpedantic)
target_include_directories(advancing_front
    PUBLIC include/NodeGenLib)
target_link_libraries(advancing_front
    Eigen3::Eigen)
target_link_libraries(advancing_front
    nanoflann::nanoflann)

# Build instructions for MEX function step_generate_boundary_nodes
matlab_add_mex(
    NAME step_generate_boundary_nodes
    SRC src/step_generate_boundary_nodes.cpp
    R2018a
)
target_compile_features(step_generate_boundary_nodes
    PUBLIC cxx_std_17)
target_compile_options(step_generate_boundary_nodes
    PUBLIC -Wall -Wextra -Wpedantic)
target_include_directories(step_generate_boundary_nodes
    PUBLIC include/NodeGenLib)
target_include_directories(step_generate_boundary_nodes
    SYSTEM PUBLIC "${OpenCASCADE_INCLUDE_DIR}")
target_link_directories(step_generate_boundary_nodes
    PUBLIC "${OpenCASCADE_LIBRARY_DIR}")
target_link_libraries(step_generate_boundary_nodes
    Eigen3::Eigen)
target_link_libraries(step_generate_boundary_nodes
    nanoflann::nanoflann)
foreach(LIB ${OpenCASCADE_LIBRARIES})
    target_link_libraries(step_generate_boundary_nodes ${LIB})
endforeach()

# Build instructions for MEX function step_plot_pcurves
matlab_add_mex(
    NAME step_plot_pcurves
    SRC src/step_plot_pcurves.cpp
    R2018a
)
target_compile_features(step_plot_pcurves
    PUBLIC cxx_std_17)
target_compile_options(step_plot_pcurves
    PUBLIC -Wall -Wextra -Wpedantic)
target_include_directories(step_plot_pcurves
    SYSTEM PUBLIC "${OpenCASCADE_INCLUDE_DIR}")
target_link_directories(step_plot_pcurves
    PUBLIC "${OpenCASCADE_LIBRARY_DIR}")
foreach(LIB ${OpenCASCADE_LIBRARIES})
    target_link_libraries(step_plot_pcurves ${LIB})
endforeach()

# Build instructions for MEX function step_print_contents
matlab_add_mex(
    NAME step_print_contents
    SRC src/step_print_contents.cpp
    R2018a
)
target_compile_features(step_print_contents
    PUBLIC cxx_std_17)
target_compile_options(step_print_contents
    PUBLIC -Wall -Wextra -Wpedantic)
target_include_directories(step_print_contents
    SYSTEM PUBLIC "${OpenCASCADE_INCLUDE_DIR}")
target_link_directories(step_print_contents
    PUBLIC "${OpenCASCADE_LIBRARY_DIR}")
foreach(LIB ${OpenCASCADE_LIBRARIES})
    target_link_libraries(step_print_contents ${LIB})
endforeach()

# Build instructions for MEX function step_print_faces_info
matlab_add_mex(
    NAME step_print_faces_info
    SRC src/step_print_faces_info.cpp
    R2018a
)
target_compile_features(step_print_faces_info
    PUBLIC cxx_std_17)
target_compile_options(step_print_faces_info
    PUBLIC -Wall -Wextra -Wpedantic)
target_include_directories(step_print_faces_info
    SYSTEM PUBLIC "${OpenCASCADE_INCLUDE_DIR}")
target_link_directories(step_print_faces_info
    PUBLIC "${OpenCASCADE_LIBRARY_DIR}")
foreach(LIB ${OpenCASCADE_LIBRARIES})
    target_link_libraries(step_print_faces_info ${LIB})
endforeach()
