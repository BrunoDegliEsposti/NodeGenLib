Installation instructions
- Install cmake
- Install the library "OpenCASCADE" for working with CAD domains
- Install the library "Eigen" for linear algebra, version 3.4 or newer
- Install the library "nanoflann" for nearest-neighbor computations
- Run "cmake ."
- Run "make"
- Add this folder to MATLAB's path

TODO:
- Test isotropic candidates around y in R^d vs isotropic candidates
around G(y) in the tangent plane spanned by the columns of dG(y)
- Remove the auto keyword from the source code
- Write a demo.cpp standalone program that discretizes a sphere
and its interior, so that performance of our advancing front method
can be analyzed with a CPU profiler
- For each accepted candidate node, store the index of its parent node
and its generation number
- Skip inclusion queries against aabb and Z if candidates belong
to a sufficiently old generation
- Use parent node of y to prune candidates around y, and additionally
use any accepted candidates to further prune candidates during expansion
