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
- Current performance is about 100k nodes/second,
and further optimizations are planned.
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

1. Install a C++ compiler, MATLAB, and the build system CMake.
2. Install the C++ library "Open CASCADE Technology",
including all possible sub-packages and development packages
For example, on Ubuntu, install all the packages with the libocct prefix.
On macOS, ```brew install opencascade``` is sufficient.
3. Install the C++ library "Eigen", version 3.4 or newer.
4. Optionally install the C++ library "nanoflann".
If nanoflann is not installed, a new copy of the library
will be automatically downloaded when building NodeGenLib.
5. Clone this repository and, using a terminal, change
the current directory to the one that you have just downloaded.
6. Run the ```cmake .``` command to generate appropriate build files.
7. Run ```make``` to build NodeGenLib.
8. Locate the newly generated MATLAB MEX functions and add
their containing folder to MATLAB's search path.

## Solution to frequent problems

- If the make command fails on a Linux distribution
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

- Write some standalone demos/tests, so that performance
of the advancing front method can be profiled and optimized.
The demo can also be used to check if NodeGenLib has been
built and installed correctly.
- For each accepted candidate node, store the index of its parent node
and its generation number. Then, skip inclusion queries against
aabb and Z if candidates belong to a sufficiently old generation.
- Use parent node of y to prune candidates around y, and additionally
use any accepted candidates to further prune candidates during expansion.
- Feel free to suggest or contribute additional ideas for further improvements
