// MEX API
//https://it.mathworks.com/help/matlab/call-mex-files-1.html
// Matrix API
//https://it.mathworks.com/help/matlab/cc-mx-matrix-library.html
// Typed Data Access
//https://it.mathworks.com/help/matlab/matlab_external/c-matrix-api-typed-data-access.html

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include "mex.h"
#include "matrix.h"
#include <STEPControl_Reader.hxx>
#include <BRepTools.hxx>
#include <BRep_Tool.hxx>
#include <Geom2d_Curve.hxx>
#include <gp_Pnt2d.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Shape.hxx>
#include <TopExp_Explorer.hxx>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Process arguments
    if (nrhs != 1) {
        mexErrMsgIdAndTxt("STEP:inputargs","Function takes in only one argument");
    }
    Standard_CString filename = mxArrayToString(prhs[0]);
    if (filename == nullptr) {
        mexErrMsgIdAndTxt("STEP:inputargs","Cannot convert filename to string");
    }

    // Read input file
    STEPControl_Reader reader;
    IFSelect_ReturnStatus read_status = reader.ReadFile(filename);
    if (read_status != IFSelect_RetDone) {
        mexErrMsgIdAndTxt("STEP:inputfile","Error reading STEP file");
    }
    Standard_Integer nShapes = reader.TransferRoots();
    mexPrintf("Input file contains %d shapes\n", nShapes);

    // Process shape
    TopoDS_Shape shape = reader.OneShape();
    switch (shape.ShapeType()) {
    case TopAbs_COMPOUND:
        mexPrintf("STEP file loaded as a COMPOUND\n");
        break;
    case TopAbs_COMPSOLID:
        mexPrintf("STEP file loaded as a COMPSOLID\n");
        break;
    case TopAbs_SOLID:
        mexPrintf("STEP file loaded as a SOLID\n");
        break;
    case TopAbs_SHELL:
        mexPrintf("STEP file loaded as a SHELL\n");
        break;
    case TopAbs_FACE:
        mexPrintf("STEP file loaded as a FACE\n");
        break;
    case TopAbs_WIRE:
        mexPrintf("STEP file loaded as a WIRE\n");
        break;
    case TopAbs_EDGE:
        mexPrintf("STEP file loaded as an EDGE\n");
        break;
    case TopAbs_VERTEX:
        mexPrintf("STEP file loaded as a VERTEX\n");
        break;
    default:
        mexErrMsgIdAndTxt("STEP:inputfile","Error loading STEP file");
    }

    // Explore shape
    TopExp_Explorer faces(shape,TopAbs_FACE);
    for (int iFace = 1; faces.More(); faces.Next(), iFace++) {
        const TopoDS_Face &face = TopoDS::Face(faces.Current());
        int nWires = 0;
        int nEdges = 0;
        Standard_Boolean allSameParameter = true;
        Standard_Boolean allSameRange = true;
        Standard_Boolean anyDegenerated = false;
        Standard_Real u0, u1, v0, v1;
        BRepTools::UVBounds(face, u0, u1, v0, v1);
        TopExp_Explorer wires(face,TopAbs_WIRE);
        for (; wires.More(); wires.Next()) {
            nWires += 1;
            const TopoDS_Wire &wire = TopoDS::Wire(wires.Current());
            TopExp_Explorer edges(wire,TopAbs_EDGE);
            for (; edges.More(); edges.Next()) {
                nEdges += 1;
                const TopoDS_Edge &edge = TopoDS::Edge(edges.Current());
                allSameParameter &= BRep_Tool::SameParameter(edge);
                allSameRange &= BRep_Tool::SameRange(edge);
                anyDegenerated |= BRep_Tool::Degenerated(edge);
            }
        }
        mexPrintf("Face %d has %d wires and %d edges\n",iFace,nWires,nEdges);
        mexPrintf("  allSameParameter: %d, allSameRange: %d, anyDegenerated: %d\n",
            int(allSameParameter),int(allSameRange),int(anyDegenerated));
        mexPrintf("  Bounding box of pcurves: (%f,%f) x (%f,%f)\n", u0,u1,v0,v1);
    }
}
