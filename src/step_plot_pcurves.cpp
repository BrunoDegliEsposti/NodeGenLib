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

void plotEdgeOnFace(int figid, const TopoDS_Face &face, int edgeid, const TopoDS_Edge &edge)
{
    // Extract pcurve
    Standard_Real t0,t1;
    Handle(Geom2d_Curve) pcurve = BRep_Tool::CurveOnSurface(edge,face,t0,t1);
    if (pcurve.IsNull()) {
        mexErrMsgIdAndTxt("STEP:nullptr","Handle to pcurve is NULL");
    }
    // Sample pcurve and its normals
    mxArray *px = mxCreateDoubleMatrix(101, 1, mxREAL);
    mxArray *py = mxCreateDoubleMatrix(101, 1, mxREAL);
    mxArray *nx = mxCreateDoubleMatrix(101, 1, mxREAL);
    mxArray *ny = mxCreateDoubleMatrix(101, 1, mxREAL);
    double *px_buf = mxGetDoubles(px);
    double *py_buf = mxGetDoubles(py);
    double *nx_buf = mxGetDoubles(nx);
    double *ny_buf = mxGetDoubles(ny);
    for (int i = 0; i <= 100; i++) {
        Standard_Real t = t0 + i*(t1-t0)/100.0;
        gp_Pnt2d p;
        gp_Vec2d v;
        pcurve->D1(t,p,v);
        px_buf[i] = p.X();
        py_buf[i] = p.Y();
        Standard_Real vNorm = v.Magnitude();
        Standard_Real s = (edge.Orientation() == face.Orientation()) ? 1.0 : -1.0;
        nx_buf[i] =  (i/100.0)*s*v.Y()/vNorm;
        ny_buf[i] = -(i/100.0)*s*v.X()/vNorm;
    }
    // Open the right figure
    mxArray *plhs_figure[0];
    mxArray *prhs_figure[1] = {mxCreateDoubleScalar(double(figid))};
    mexCallMATLAB(0, plhs_figure, 1, prhs_figure, "figure");
    // Plot pcurve points
    const char* colors[7] = {"r", "g", "b", "c", "m", "y", "k"};
    mxArray *plhs_plot[0];
    mxArray *prhs_plot[3] = {px,py,mxCreateString(colors[(edgeid-1)%7])};
    mexCallMATLAB(0, plhs_plot, 3, prhs_plot, "plot");
    // Call hold on
    mxArray *plhs_hold[0];
    mxArray *prhs_hold[1] = {mxCreateString("on")};
    mexCallMATLAB(0, plhs_hold, 1, prhs_hold, "hold");
    // Plot pcurve normals
    mxArray *plhs_quiver[0];
    mxArray *prhs_quiver[4] = {px,py,nx,ny};
    mexCallMATLAB(0, plhs_quiver, 4, prhs_quiver, "quiver");
}

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
    if (nShapes == 0) {
        mexErrMsgIdAndTxt("STEP:inputfile","Error loading STEP file");
    }
    TopoDS_Shape shape = reader.OneShape();

    const char *orientation_names[4] = {
        "TopAbs_FORWARD", "TopAbs_REVERSED", "TopAbs_INTERNAL", "TopAbs_EXTERNAL"
    };

    // Explore shape
    TopExp_Explorer faces(shape,TopAbs_FACE);
    for (int i = 0; faces.More(); faces.Next(), i++) {
        const TopoDS_Face &face = TopoDS::Face(faces.Current());
        mexPrintf("Face %d found\n",i+1);
        mexPrintf("Face has orientation %s\n",orientation_names[face.Orientation()]);
        TopExp_Explorer wires(face,TopAbs_WIRE);
        for (int j = 0; wires.More(); wires.Next(), j++) {
            const TopoDS_Wire &wire = TopoDS::Wire(wires.Current());
            mexPrintf("  Wire %d found on face %d\n",j+1,i+1);
            mexPrintf("  Wire has orientation %s\n",orientation_names[wire.Orientation()]);
            TopExp_Explorer edges(wire,TopAbs_EDGE);
            for (int k = 0; edges.More(); edges.Next(), k++) {
                const TopoDS_Edge &edge = TopoDS::Edge(edges.Current());
                mexPrintf("    Edge %d found on wire %d on face %d\n",k+1,j+1,i+1);
                mexPrintf("    Edge has orientation %s\n",orientation_names[edge.Orientation()]);
                plotEdgeOnFace(i+1, face, k+1, edge);
            }
        }
    }
}
