/*
 * NodeGenLib - A library for advancing front node generation
 * Copyright (C) 2024 Bruno Degli Esposti
 *
 * This file is part of NodeGenLib.
 *
 * NodeGenLib is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * NodeGenLib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with NodeGenLib. If not, see <http://www.gnu.org/licenses/>.
 */

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
#include <algorithm>
#include "mex.h"
#include "matrix.h"
#include <Eigen/Dense>
#include <nanoflann.hpp>
#include "advancing_front.hpp"

#include <STEPControl_Reader.hxx>
#include <BRepTools.hxx>
#include <BRep_Tool.hxx>
#include <Geom2d_Curve.hxx>
#include <Geom2dAdaptor_Curve.hxx>
#include <gp_Pnt2d.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Shape.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>

void plotParametricNodes(const Nodes<2> &Z, const Nodes<2> &Y)
{
    // Sample pcurve and its normals
    size_t NZ = Z.points.size();
    size_t NY = Y.points.size();
    mxArray *Zpx_mxArray = mxCreateDoubleMatrix(NZ, 1, mxREAL);
    mxArray *Zpy_mxArray = mxCreateDoubleMatrix(NZ, 1, mxREAL);
    mxArray *Znx_mxArray = mxCreateDoubleMatrix(NZ, 1, mxREAL);
    mxArray *Zny_mxArray = mxCreateDoubleMatrix(NZ, 1, mxREAL);
    mxArray *Ypx_mxArray = mxCreateDoubleMatrix(NY, 1, mxREAL);
    mxArray *Ypy_mxArray = mxCreateDoubleMatrix(NY, 1, mxREAL);
    double *Zpx_buf = mxGetDoubles(Zpx_mxArray);
    double *Zpy_buf = mxGetDoubles(Zpy_mxArray);
    double *Znx_buf = mxGetDoubles(Znx_mxArray);
    double *Zny_buf = mxGetDoubles(Zny_mxArray);
    double *Ypx_buf = mxGetDoubles(Ypx_mxArray);
    double *Ypy_buf = mxGetDoubles(Ypy_mxArray);
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
    // Plot Z.points
    mxArray *ro_mxArray = mxCreateString("ro");
    mxArray *prhs_Zplot[3] = {Zpx_mxArray,Zpy_mxArray,ro_mxArray};
    mexCallMATLAB(0, nullptr, 3, prhs_Zplot, "plot");
    mxDestroyArray(ro_mxArray);
    // Call hold on
    mxArray *on_mxArray = mxCreateString("on");
    mxArray *prhs_holdon[1] = {on_mxArray};
    mexCallMATLAB(0, nullptr, 1, prhs_holdon, "hold");
    mxDestroyArray(on_mxArray);
    // Plot Z.normals
    mxArray *prhs_quiver[4] = {Zpx_mxArray,Zpy_mxArray,Znx_mxArray,Zny_mxArray};
    mexCallMATLAB(0, nullptr, 4, prhs_quiver, "quiver");
    // Plot Y.points
    mxArray *kdot_mxArray = mxCreateString("k.");
    mxArray *prhs_Yplot[3] = {Ypx_mxArray,Ypy_mxArray,kdot_mxArray};
    mexCallMATLAB(0, nullptr, 3, prhs_Yplot, "plot");
    mxDestroyArray(kdot_mxArray);
    // Call hold off
    mxArray *off_mxArray = mxCreateString("off");
    mxArray *prhs_holdoff[1] = {off_mxArray};
    mexCallMATLAB(0, nullptr, 1, prhs_holdoff, "hold");
    mxDestroyArray(off_mxArray);
    // Call axis equal
    mxArray *equal_mxArray = mxCreateString("equal");
    mxArray *prhs_axis[1] = {equal_mxArray};
    mexCallMATLAB(0, nullptr, 1, prhs_axis, "axis");
    mxDestroyArray(equal_mxArray);
    // Free up memory
    mxDestroyArray(Zpx_mxArray);
    mxDestroyArray(Zpy_mxArray);
    mxDestroyArray(Znx_mxArray);
    mxDestroyArray(Zny_mxArray);
    mxDestroyArray(Ypx_mxArray);
    mxDestroyArray(Ypy_mxArray);
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
// function [points,normals] = step_generate_boundary_nodes(filename,iSolid,h,Nmax,seed,iFaceDebug)
{
    // Rename input and output arguments
    if (nrhs != 6) {
        mexErrMsgIdAndTxt("StepGenerateBoundaryNodes:Input","Wrong number of input arguments");
    }
    if (nlhs != 2) {
        mexErrMsgIdAndTxt("StepGenerateBoundaryNodes:Input","Wrong number of output arguments");
    }
    const mxArray *filename_mxArray = prhs[0];
    const mxArray *iSolid_mxArray = prhs[1];
    const mxArray *h_mxArray = prhs[2];
    const mxArray *Nmax_mxArray = prhs[3];
    const mxArray *seed_mxArray = prhs[4];
    const mxArray *iFaceDebug_mxArray = prhs[5];
    mxArray *points_mxArray = nullptr;
    mxArray *normals_mxArray = nullptr;
    
    // Unwrap filename
    Standard_CString filename = mxArrayToString(filename_mxArray);
    if (filename == nullptr) {
        mexErrMsgIdAndTxt("StepGenerateBoundaryNodes:Input","Filename must be a char array");
    }
    
    // Unwrap iSolid
    if (!mxIsScalar(iSolid_mxArray) || !mxIsNumeric(iSolid_mxArray)) {
        mexErrMsgIdAndTxt("StepGenerateBoundaryNodes:Input","iSolid must be a scalar of numeric type");
    }
    int iSolid = (int)mxGetScalar(iSolid_mxArray);
    
    // Unwrap h
    if (mxGetClassID(h_mxArray) != mxFUNCTION_CLASS) {
        mexErrMsgIdAndTxt("StepGenerateBoundaryNodes:Input","h is not a function handle");
    }
    auto h = [h_mxArray](const Point<3> &x) -> double {
        mxArray *plhs_feval[1];
        mxArray *prhs_feval[2];
        prhs_feval[0] = (mxArray*)h_mxArray;
        prhs_feval[1] = mxArray_from_point<3>(x);
        int errcode = mexCallMATLAB(1, plhs_feval, 2, prhs_feval, "feval");
        if (errcode) {
            mexErrMsgIdAndTxt("StepGenerateBoundaryNodes:Feval",
                              "Evaluation of h produced errors");
        }
        if (!mxIsDouble(plhs_feval[0]) || mxGetNumberOfElements(plhs_feval[0]) != 1) {
            mexErrMsgIdAndTxt("StepGenerateBoundaryNodes:Feval",
                              "Output of h must be a scalar");
        }
        double hx = mxGetScalar(plhs_feval[0]);
        mxDestroyArray(plhs_feval[0]);
        mxDestroyArray(prhs_feval[1]);
        return hx;
    };
    
    // Unwrap Nmax
    if (!mxIsScalar(Nmax_mxArray) || !mxIsNumeric(Nmax_mxArray)) {
        mexErrMsgIdAndTxt("StepGenerateBoundaryNodes:Input","Nmax must be a scalar of numeric type");
    }
    size_t Nmax = (size_t)mxGetScalar(Nmax_mxArray);
    size_t Z_size_limit = std::max((size_t)1000000,Nmax);

    // Unwrap seed
    if (!mxIsScalar(seed_mxArray) || !mxIsNumeric(seed_mxArray)) {
        mexErrMsgIdAndTxt("StepGenerateBoundaryNodes:Input","seed must be a scalar of numeric type");
    }
    uint64_t seed = (uint64_t)mxGetScalar(seed_mxArray);
    std::mt19937_64 rng(seed);
    
    // Unwrap iFaceDebug
    if (!mxIsScalar(iFaceDebug_mxArray) || !mxIsNumeric(iFaceDebug_mxArray)) {
        mexErrMsgIdAndTxt("StepGenerateBoundaryNodes:Input","iFaceDebug must be a scalar of numeric type");
    }
    int iFaceDebug = (int)mxGetScalar(iFaceDebug_mxArray);

    // Read input file
    STEPControl_Reader reader;
    IFSelect_ReturnStatus read_status = reader.ReadFile(filename);
    if (read_status != IFSelect_RetDone) {
        mexErrMsgIdAndTxt("StepGenerateBoundaryNodes:Input","Error reading STEP file");
    }
    Standard_Integer nShapes = reader.TransferRoots();
    if (nShapes == 0) {
        mexErrMsgIdAndTxt("StepGenerateBoundaryNodes:Input","Error loading STEP file");
    }
    TopoDS_Shape shape = reader.OneShape();
    
    // Focus on the solid with index iSolid in the loaded shape.
    // Indices in TopTools_IndexedMapOfShape start from 1, not from 0
    TopTools_IndexedMapOfShape solids_map;
    TopExp::MapShapes(shape, TopAbs_SOLID, solids_map);
    Standard_Integer nSolids = solids_map.Extent();
    if (iSolid < 1 || iSolid > nSolids) {
        mexErrMsgIdAndTxt("StepGenerateBoundaryNodes:Input","Solid index out of range [1,%d]",nSolids);
    }

    // Variable GY stores the generated boundary nodes on the solid.
    // Grow GY by iterating over all faces in the solid
    Nodes<3> GY;
    TopExp_Explorer faces(solids_map(iSolid),TopAbs_FACE);
    for (int i = 0; faces.More() && Nmax > 0; faces.Next(), i++)
    {
        // Fetch parametrization F and bounding box of parametric domain
        const TopoDS_Face &face = TopoDS::Face(faces.Current());
        Handle(Geom_Surface) F_handle = BRep_Tool::Surface(face);
        double U1, U2, V1, V2;
        BRepTools::UVBounds(face,U1,U2,V1,V2);
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
        // Grow GY by iterating over all edges in this face and update hmin
        double hmin = std::numeric_limits<double>::infinity();
        TopExp_Explorer edges(face,TopAbs_EDGE);
        for (int j = 0; edges.More() && Nmax > 0; edges.Next(), j++)
        {
            // Fetch parametrization g: [t1,t2] -> R^2
            const TopoDS_Edge &edge = TopoDS::Edge(edges.Current());
            double t1, t2;
            Handle(Geom2d_Curve) g_handle = BRep_Tool::CurveOnSurface(edge,face,t1,t2);
            Geom2dAdaptor_Curve g_adaptor(g_handle,t1,t2);
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
            // Discretize the interval [t1,t2] according to h: R^3 -> R
            AABB<1> aabb_G(t1,t2);
            Nodes<1> t_ij;
            advancing_front<1,3>(aabb_G, Nodes<1>(), t_ij, G, dG, h, Nmax, rng);
            size_t Nt = t_ij.points.size();
            // Discretize the image of g in parametric domain
            for (size_t k = 0; k < Nt; k++) {
                Point<1> t = t_ij.points[k];
                Point<2> gt = g(t);
                double s_edge = (edge.Orientation() == face.Orientation()) ? 1.0 : -1.0;
                Point<2> nu_edge = normal_from_jacobian<1,2>(dg(t)) * s_edge;
                double etol = BRep_Tool::Tolerance(edge);
                double alpha = 2.5 * etol / ((dF(gt)*nu_edge).norm() + 1e-15);
                Point<3> Fnew = F(gt - alpha*nu_edge);
                Eigen::Matrix<double,3,2> dFnew = dF(gt - alpha*nu_edge);
                GY.points.push_back(Fnew);
                double s_face = face.Orientation() ? -1.0 : 1.0;
                GY.normals.push_back(normal_from_jacobian(dFnew) * s_face);
                Nmax = Nmax - 1;
            }
            // Find hmin as the minimum distance between consecutive nodes in parameter space
            for (size_t k = 0; k < Nt-1; k++) {
                double dist = (g(t_ij.points[k+1])-g(t_ij.points[k])).norm();
                hmin = std::min(hmin,dist);
            }
        }
        // If hmin has not been updated, skip further processing because the face is too small
        if (hmin == std::numeric_limits<double>::infinity()) {
            continue;
        }
        // Grow Z_i by iterating over all edges in this face
        Nodes<2> Z_i;
        edges.ReInit();
        for (int j = 0; edges.More() && Nmax > 0; edges.Next(), j++)
        {
            // Fetch parametrization g: [t1,t2] -> R^2
            const TopoDS_Edge &edge = TopoDS::Edge(edges.Current());
            double t1, t2;
            Handle(Geom2d_Curve) g_handle = BRep_Tool::CurveOnSurface(edge,face,t1,t2);
            Geom2dAdaptor_Curve g_adaptor(g_handle,t1,t2);
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
            // Discretize the interval [t1,t2] according to constant spacing hmin.
            // If the edge is too short, skip the discretization process
            AABB<1> aabb_g(t1,t2);
            Nodes<1> t_ij;
            Point<1> tL(t1 + hmin/(2*dg(Point<1>(t1)).norm()+1e-15));
            Point<1> tR(t2 - hmin/(2*dg(Point<1>(t2)).norm()+1e-15));
            if (tL(0) >= tR(0)) {
                continue;
            }
            t_ij.points.push_back(tL);
            t_ij.points.push_back(tR);
            Nodes<2> gt = advancing_front<1,2>(aabb_g, Nodes<1>(), t_ij, g, dg,
                                               [hmin](const Point<2>) -> double {return hmin;},
                                               Z_size_limit, rng);
            size_t Ngt = gt.points.size();
            if (Ngt == Z_size_limit) {
                mexErrMsgIdAndTxt("StepGenerateBoundaryNodes:Discretization",
                  "Number of nodes on trimming curve %d on face %d exceeded Z_size_limit", j+1, i+1);
            }
            double s = (edge.Orientation() == face.Orientation()) ? 1.0 : -1.0;
            for (size_t k = 0; k < Ngt; k++) {
                Z_i.points.push_back(gt.points[k]);
                Z_i.normals.push_back(gt.normals[k] * s);
            }
            Point<1> t_last1 = t_ij.points[Ngt-1];
            Point<1> t_last2 = t_ij.points[Ngt-2];
            Point<2> gt_last1 = gt.points[Ngt-1];
            Point<2> gt_last2 = gt.points[Ngt-2];
            if ((gt_last1-gt_last2).norm() > 1.2*hmin) {
                Point<1> t_last = 0.5*t_last1 + 0.5*t_last2;
                Normal<2> nu_last = normal_from_jacobian(dg(t_last));
                Z_i.points.push_back(g(t_last));
                Z_i.normals.push_back(nu_last * s);
            }
        }
        // Use the uniform boundary nodes Z_i to place points on the face
        Nodes<2> Y_i;
        Nodes<3> GY_i = advancing_front<2,3>(aabb_F, Z_i, Y_i, F, dF, h, Nmax + Z_i.points.size(), rng);
        // Append GY_i to GY skipping the first |Z_i| elements
        for (size_t k = Z_i.points.size(); k < GY_i.points.size(); k++) {
            GY.points.push_back(GY_i.points[k]);
            double s = face.Orientation() ? -1.0 : 1.0;
            GY.normals.push_back(GY_i.normals[k] * s);
            Nmax = Nmax - 1;
        }
        // Plot parametric nodes if debug is enabled
        if (i+1 == iFaceDebug) {
            plotParametricNodes(Z_i,Y_i);
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
