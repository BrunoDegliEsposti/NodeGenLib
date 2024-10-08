# NodeGenLib

**NodeGenLib** is a header-only C++ library
and a collection of MATLAB MEX functions
to generate nodes on complex 2D and 3D
domains using an advancing front method. 
Its aim is to make it easier to run numerical
simulations and solve practical problems
using meshless/meshfree methods
on challenging geometries.

## Features

- Nodes can be generated on any kind of domain,
as long as a parametrization of its boundary is known,
possibly consisting of multiple patches.
Nodes are generated on each patch first, and then
more nodes are generated in the interior.
- Integration with Open CASCADE Technology allows
3D domains to be given in B-Rep format using STEP files.
Trimmed geometries are supported.
- Current performance is about 5k nodes/second in 3D using MATLAB.
Easy opportunities for optimization have not been
taken advantage of yet, although this is planned for the near future.
- Variable node density is supported. The local node
spacing around each point is prescribed by a scalar field
passed as a lambda expression (C++) or function handle (MATLAB).
- To prevent nodes from being generated outside of the
domain, a fast and meshless point membership classification
algorithm is used. In the case of CAD geometries, the same
algorithm is also used to prevent nodes from being generated
outside of the trimmed parameter domain of each patch.

## Installation

To install NodeGenLib, follow these steps:

1. Install a compiler that supports C++17.
2. Install the build system CMake, minimum version 3.12.
3. Install MATLAB, minimum version R2018a.
4. Install the C++ library "Open CASCADE Technology",
minimum version 7.7.0. Include all possible sub-packages
and development packages. For example, on Ubuntu,
install all the packages with the libocct prefix.
On macOS, ```brew install opencascade``` is enough.
5. Install the C++ library "Eigen", minimum version 3.4.
6. Clone this repository and, using a terminal, change
the current directory to it.
7. Create a build directory inside NodeGenLib
and move there: ```mkdir build && cd build```.
8. Execute the command ```cmake -D CMAKE_BUILD_TYPE=Release ..```
to generate appropriate build files.
9. Run ```cmake --build .``` to build NodeGenLib.
10. Locate the newly generated MEX functions inside the build
folder (or a subfolder called Release). Add the folder containing
the MEX functions to MATLAB's search path.

## Solution to frequent problems

- If the build command fails on a Linux distribution
with an error related to fontconfig,
install the development package for fontconfig,
such as "fontconfig-devel".

## License

This project is licensed under version 3 of the LGPL license.
See the LICENSE file for details.

## Contributing

Pull requests are welcome, but please contact the main author
first at brunodegliesposti@gmail.com to make sure that your
contribution falls within the scope of the library and
fits its design guidelines.

## Ideas for future work

- Write documentation to describe how NodeGenLib can be used
as a header-only C++ library.
- Write documentation to describe the MATLAB MEX functions
provided by NodeGenLib, their input arguments and their
output arguments.
- Write standalone demos and tests, so that performance
of the advancing front method can be profiled and optimized.
Tests are also needed to check if NodeGenLib has been
built and installed correctly on a given system.
- For each accepted candidate node, store the index of its parent node
and its generation number (number of parents).
Then, skip inclusion queries if candidates belong to a sufficiently old generation.
- Use parent of node y to prune candidates around y, and additionally
use any accepted candidates to further prune candidates during expansion.
- Feel free to suggest or contribute additional ideas for further improvements
