// MEX API
//https://it.mathworks.com/help/matlab/call-mex-files-1.html
// Matrix API
//https://it.mathworks.com/help/matlab/cc-mx-matrix-library.html
// Typed Data Access
//https://it.mathworks.com/help/matlab/matlab_external/c-matrix-api-typed-data-access.html

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <cmath>
#include <random>
#include "mex.h"
#include "matrix.h"
#include <Eigen/Dense>
#include <nanoflann.hpp>
#include "advancing_front.hpp"

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

void plotParametricNodes(int figid, const Nodes<2> &Z, const Nodes<2> &Y)
{
    // Sample pcurve and its normals
    size_t NZ = Z.points.size();
    size_t NY = Y.points.size();
    mxArray *Zpx = mxCreateDoubleMatrix(NZ, 1, mxREAL);
    mxArray *Zpy = mxCreateDoubleMatrix(NZ, 1, mxREAL);
    mxArray *Znx = mxCreateDoubleMatrix(NZ, 1, mxREAL);
    mxArray *Zny = mxCreateDoubleMatrix(NZ, 1, mxREAL);
    mxArray *Ypx = mxCreateDoubleMatrix(NY, 1, mxREAL);
    mxArray *Ypy = mxCreateDoubleMatrix(NY, 1, mxREAL);
    double *Zpx_buf = mxGetDoubles(Zpx);
    double *Zpy_buf = mxGetDoubles(Zpy);
    double *Znx_buf = mxGetDoubles(Znx);
    double *Zny_buf = mxGetDoubles(Zny);
    double *Ypx_buf = mxGetDoubles(Ypx);
    double *Ypy_buf = mxGetDoubles(Ypy);
    for (size_t i = 0; i <= NZ; i++) {
        Zpx_buf[i] = Z.points[i][0];
        Zpy_buf[i] = Z.points[i][1];
        Znx_buf[i] = Z.normals[i][0];
        Zny_buf[i] = Z.normals[i][1];
    }
    for (size_t i = 0; i <= NY; i++) {
        Ypx_buf[i] = Y.points[i][0];
        Ypy_buf[i] = Y.points[i][1];
    }
    // Open the right figure
    mxArray *prhs_figure[1] = {mxCreateDoubleScalar(double(figid))};
    mexCallMATLAB(0, nullptr, 1, prhs_figure, "figure");
    // Plot Z.points
    mxArray *prhs_Zplot[3] = {Zpx,Zpy,mxCreateString("ro")};
    mexCallMATLAB(0, nullptr, 3, prhs_Zplot, "plot");
    // Call hold on
    mxArray *prhs_holdon[1] = {mxCreateString("on")};
    mexCallMATLAB(0, nullptr, 1, prhs_holdon, "hold");
    // Plot Z.normals
    mxArray *prhs_quiver[4] = {Zpx,Zpy,Znx,Zny};
    mexCallMATLAB(0, nullptr, 4, prhs_quiver, "quiver");
    // Plot Y.points
    mxArray *prhs_Yplot[3] = {Ypx,Ypy,mxCreateString("k.")};
    mexCallMATLAB(0, nullptr, 3, prhs_Yplot, "plot");
    // Call hold off
    mxArray *prhs_holdoff[1] = {mxCreateString("off")};
    mexCallMATLAB(0, nullptr, 1, prhs_holdoff, "hold");
    // Call axis equal
    mxArray *prhs_axis[1] = {mxCreateString("equal")};
    mexCallMATLAB(0, nullptr, 1, prhs_axis, "axis");
    // Free up memory
    mxDestroyArray(Zpx);
    mxDestroyArray(Zpy);
    mxDestroyArray(Znx);
    mxDestroyArray(Zny);
    mxDestroyArray(Ypx);
    mxDestroyArray(Ypy);
}

template <int dim>
mxArray* mxArray_from_point(const Point<dim> &p)
{
    mxArray *p_mxArray = mxCreateDoubleMatrix(1,dim,mxREAL);
    double *p_buf = mxGetDoubles(p_mxArray);
    for (size_t i = 0; i < dim; i++) {
        p_buf[i] = p(i);
    }
    return p_mxArray;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
// function [points,normals] = step_advancing_front(filename,h,Nmax,seed,debug)
{
    // Rename input and output arguments
    if (nrhs != 5) {
        mexErrMsgIdAndTxt("StepAdvancingFront:Input","Wrong number of input arguments");
    }
    if (nlhs != 2) {
        mexErrMsgIdAndTxt("StepAdvancingFront:Input","Wrong number of output arguments");
    }
    const mxArray *filename_mxArray = prhs[0];
    const mxArray *h_mxArray = prhs[1];
    const mxArray *Nmax_mxArray = prhs[2];
    const mxArray *seed_mxArray = prhs[3];
    const mxArray *debug_mxArray = prhs[4];
    mxArray *points_mxArray = nullptr;
    mxArray *normals_mxArray = nullptr;
    
    // Unwrap filename
    Standard_CString filename = mxArrayToString(filename_mxArray);
    if (filename == nullptr) {
        mexErrMsgIdAndTxt("StepAdvancingFront:Input","Filename must be a char array");
    }
    
    // Unwrap h
    if (mxGetClassID(h_mxArray) != mxFUNCTION_CLASS) {
        mexErrMsgIdAndTxt("StepAdvancingFront:Input","h is not a function handle");
    }
    auto h = [h_mxArray](const Point<3> &x) -> double {
        mxArray *plhs_feval[1];
        mxArray *prhs_feval[2];
        prhs_feval[0] = (mxArray*)h_mxArray;
        prhs_feval[1] = mxArray_from_point<3>(x);
        int errcode = mexCallMATLAB(1, plhs_feval, 2, prhs_feval, "feval");
        if (errcode) {
            mexErrMsgIdAndTxt("StepAdvancingFront:Feval",
                              "Evaluation of h produced errors");
        }
        if (!mxIsDouble(plhs_feval[0]) || mxGetNumberOfElements(plhs_feval[0]) != 1) {
            mexErrMsgIdAndTxt("StepAdvancingFront:Feval",
                              "Output of h must be a scalar");
        }
        double hx = mxGetScalar(plhs_feval[0]);
        mxDestroyArray(plhs_feval[0]);
        mxDestroyArray(prhs_feval[1]);
        return hx;
    };
    
    // Unwrap Nmax
    if (!mxIsScalar(Nmax_mxArray) || !mxIsNumeric(Nmax_mxArray)) {
        mexErrMsgIdAndTxt("StepAdvancingFront:Input","Nmax must be a scalar of numeric type");
    }
    size_t Nmax = (size_t)mxGetScalar(Nmax_mxArray);

    // Unwrap seed
    if (!mxIsScalar(seed_mxArray) || !mxIsNumeric(seed_mxArray)) {
        mexErrMsgIdAndTxt("StepAdvancingFront:Input","seed must be a scalar of numeric type");
    }
    uint64_t seed = (uint64_t)mxGetScalar(seed_mxArray);
    std::mt19937_64 rng(seed);
    
    // Unwrap debug
    if (!mxIsScalar(debug_mxArray) || !mxIsLogical(debug_mxArray)) {
        mexErrMsgIdAndTxt("StepAdvancingFront:Input","debug must be a scalar of logical type");
    }
    bool debug = mxIsLogicalScalarTrue(debug_mxArray);

    // Read input file
    STEPControl_Reader reader;
    IFSelect_ReturnStatus read_status = reader.ReadFile(filename);
    if (read_status != IFSelect_RetDone) {
        mexErrMsgIdAndTxt("StepAdvancingFront:Input","Error reading STEP file");
    }
    Standard_Integer nShapes = reader.TransferRoots();
    if (nShapes == 0) {
        mexErrMsgIdAndTxt("StepAdvancingFront:Input","Error loading STEP file");
    }
    TopoDS_Shape shape = reader.OneShape();

    // Grow GY by iterating over all faces in the shape
    Nodes<3> GY;
    TopExp_Explorer faces(shape,TopAbs_FACE);
    for (int i = 0; faces.More(); faces.Next(), i++)
    {
        // Fetch parametrization F and bounding box of parametric domain
        const TopoDS_Face &face = TopoDS::Face(faces.Current());
        Handle(Geom_Surface) F_handle = BRep_Tool::Surface(face);
        double U1, U2, V1, V2;
        F_handle->Bounds(U1,U2,V1,V2);
        Point<2> blc = {U1,V1};
        Point<2> trc = {U2,V2};
        AABB<2> aabb_F(blc,trc);
        auto F = [F_handle](const Point<2> &p) -> Point<3> {
            gp_Pnt Fp;
            F_handle->D0(p(0), p(1), Fp);
            return {Fp.X(), Fp.Y(), Fp.Z()};
        };
        auto dF = [F_handle](const Point<2> &p) -> Eigen::Matrix<double,3,2> {
            gp_Pnt Fp;
            gp_Vec dFp1, dFp2;
            F_handle->D1(p(0), p(1), Fp, dFp1, dFp2);
            Eigen::Matrix<double,3,2> J {
                {dFp1.X(), dFp2.X()},
                {dFp1.Y(), dFp2.Y()},
                {dFp1.Z(), dFp2.Z()}};
            return J;
        };
        // Grow Z_i by iterating over all edges in this face
        Nodes<2> Z_i;
        TopExp_Explorer edges(face,TopAbs_EDGE);
        for (int j = 0; edges.More(); edges.Next(), j++)
        {
            // Fetch parametrization g: [t1,t2] -> R^2
            const TopoDS_Edge &edge = TopoDS::Edge(edges.Current());
            double t1, t2;
            Handle(Geom2d_Curve) g_handle = BRep_Tool::CurveOnSurface(edge,face,t1,t2);
            auto g = [g_handle](const Point<1> &t) -> Point<2> {
                gp_Pnt2d gt;
                g_handle->D0(t(0),gt);
                return {gt.X(), gt.Y()};
            };
            auto dg = [g_handle](const Point<1> &t) -> Point<2> {
                gp_Pnt2d gt;
                gp_Vec2d dgt;
                g_handle->D1(t(0),gt,dgt);
                return {dgt.X(), dgt.Y()};
            };
            // Assemble parametrization G: [t1,t2] -> R^3 by composition
            auto G = [F,g](const Point<1> &t) -> Point<3> {
                return F(g(t));
            };
            auto dG = [dF,g,dg](const Point<1> &t) -> Point<3> {
                return dF(g(t)) * dg(t);
            };
            // Discretize the interval [t1,t2] according to spacing h
            AABB<1> aabb_G(t1,t2);
            Nodes<1> t_ij;
            t_ij.points.push_back(Point<1>(t1+1e-8*(t2-t1)));
            t_ij.points.push_back(Point<1>(t2-1e-8*(t2-t1)));
            advancing_front<1,3>(aabb_G, Nodes<1>(), t_ij, G, dG, h, Nmax, rng);
            size_t Nt = t_ij.points.size();
            // Discretize the image of g in parametric domain
            for (size_t k = 0; k < Nt; k++) {
                Point<1> t = t_ij.points[k];
                double s = (edge.Orientation() == face.Orientation()) ? 1.0 : -1.0;
                Z_i.points.push_back(g(t));
                Z_i.normals.push_back(normal_from_jacobian<1,2>(dg(t)) * s);
            }
        }
        // Use Z_i to place points on the face
        Nodes<2> Y_i;
        Nodes<3> GY_i = advancing_front<2,3>(aabb_F, Z_i, Y_i, F, dF, h, Nmax, rng);
        // Append GY_i to GY skipping the first |Z_i| elements
        for (size_t k = Z_i.points.size(); k < GY_i.points.size(); k++) {
            GY.points.push_back(GY_i.points[k]);
            double s = face.Orientation() ? -1.0 : 1.0;
            GY.normals.push_back(GY_i.normals[k] * s);
        }
        // Plot parametric nodes if debug is enabled
        if (debug) {
            plotParametricNodes(i+1,Z_i,Y_i);
        }
    }
    
    // Output GY.points
    size_t Npoints = GY.points.size();
    points_mxArray = mxCreateDoubleMatrix(Npoints, 3, mxREAL);
    double *points_buf = mxGetDoubles(points_mxArray);
    for (size_t i = 0; i < Npoints; i++) {
        for (size_t j = 0; j < 3; j++) {
            points_buf[i+j*Npoints] = GY.points[i][j];
        }
    }
    
    // Output GY.normals
    size_t Nnormals = GY.normals.size();
    normals_mxArray = mxCreateDoubleMatrix(Nnormals, 3, mxREAL);
    double *normals_buf = mxGetDoubles(normals_mxArray);
    for (size_t i = 0; i < Nnormals; i++) {
        for (size_t j = 0; j < 3; j++) {
            normals_buf[i+j*Nnormals] = GY.normals[i][j];
        }
    }
    
    plhs[0] = points_mxArray;
    plhs[1] = normals_mxArray;
}
