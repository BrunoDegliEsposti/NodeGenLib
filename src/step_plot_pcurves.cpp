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
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>

void plotEdgeOnFace(const TopoDS_Face &face, int iEdge, const TopoDS_Edge &edge)
{
    // Extract pcurve
    Standard_Real t0,t1;
    Handle(Geom2d_Curve) pcurve = BRep_Tool::CurveOnSurface(edge,face,t0,t1);
    if (pcurve.IsNull()) {
        mexErrMsgIdAndTxt("StepPlotPcurves:nullptr","Handle to pcurve is NULL");
    }
    // Sample pcurve and its normals
    mxArray *px = mxCreateDoubleMatrix(101, 1, mxREAL);
    mxArray *py = mxCreateDoubleMatrix(101, 1, mxREAL);
    mxArray *qx = mxCreateDoubleMatrix( 11, 1, mxREAL);
    mxArray *qy = mxCreateDoubleMatrix( 11, 1, mxREAL);
    mxArray *nx = mxCreateDoubleMatrix( 11, 1, mxREAL);
    mxArray *ny = mxCreateDoubleMatrix( 11, 1, mxREAL);
    double *px_buf = mxGetDoubles(px);
    double *py_buf = mxGetDoubles(py);
    double *qx_buf = mxGetDoubles(qx);
    double *qy_buf = mxGetDoubles(qy);
    double *nx_buf = mxGetDoubles(nx);
    double *ny_buf = mxGetDoubles(ny);
    for (int i = 0; i <= 100; i++) {
        Standard_Real t = t0 + i*(t1-t0)/100.0;
        gp_Pnt2d p = pcurve->Value(t);
        px_buf[i] = p.X();
        py_buf[i] = p.Y();
    }
    for (int i = 0; i <= 10; i++) {
        Standard_Real t = t0 + i*(t1-t0)/10.0;
        gp_Pnt2d q;
        gp_Vec2d v;
        pcurve->D1(t,q,v);
        qx_buf[i] = q.X();
        qy_buf[i] = q.Y();
        Standard_Real vNorm = v.Magnitude();
        Standard_Real s = (edge.Orientation() == face.Orientation()) ? 1.0 : -1.0;
        nx_buf[i] =  s*v.Y()/vNorm;
        ny_buf[i] = -s*v.X()/vNorm;
    }
    // Plot pcurve points
    const char* colors[2] = {"r", "b"};
    mxArray *color_mxArray = mxCreateString(colors[iEdge%2]);
    mxArray *prhs_plot[3] = {px,py,color_mxArray};
    mexCallMATLAB(0, nullptr, 3, prhs_plot, "plot");
    // Call hold on
    mxArray *prhs_hold[1] = {mxCreateString("on")};
    mexCallMATLAB(0, nullptr, 1, prhs_hold, "hold");
    // Plot pcurve normals
    mxArray *prhs_quiver[5] = {qx,qy,nx,ny,color_mxArray};
    mexCallMATLAB(0, nullptr, 5, prhs_quiver, "quiver");
    // Free all the allocated memory
    mxDestroyArray(color_mxArray);
    mxDestroyArray(px);
    mxDestroyArray(py);
    mxDestroyArray(qx);
    mxDestroyArray(qy);
    mxDestroyArray(nx);
    mxDestroyArray(ny);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
// function [] = step_plot_pcurves(filename,iFace)
{
    // Rename input arguments
    if (nrhs != 2) {
        mexErrMsgIdAndTxt("StepPlotPcurves:Input","Function takes in two arguments");
    }
    const mxArray *filename_mxArray = prhs[0];
    const mxArray *iFace_mxArray = prhs[1];
    
    // Unwrap filename
    Standard_CString filename = mxArrayToString(filename_mxArray);
    if (filename == nullptr) {
        mexErrMsgIdAndTxt("StepPlotPcurves:Input","Filename must be a char array");
    }
    
    // Unwrap iFace
    if (!mxIsScalar(iFace_mxArray) || !mxIsNumeric(iFace_mxArray)) {
        mexErrMsgIdAndTxt("StepPlotPcurves:Input","Face index must be a scalar of numeric type");
    }
    Standard_Integer iFace = (Standard_Integer)mxGetScalar(iFace_mxArray);

    // Read input file
    STEPControl_Reader reader;
    IFSelect_ReturnStatus read_status = reader.ReadFile(filename);
    if (read_status != IFSelect_RetDone) {
        mexErrMsgIdAndTxt("StepPlotPcurves:Input","Error reading STEP file");
    }
    Standard_Integer nShapes = reader.TransferRoots();
    if (nShapes == 0) {
        mexErrMsgIdAndTxt("StepPlotPcurves:Input","Error loading STEP file");
    }
    TopoDS_Shape shape = reader.OneShape();
    
    // Make a list of all faces, with no duplicates.
    // Indices of TopTools_IndexedMapOfShape start from 1, not from 0.
    TopTools_IndexedMapOfShape faces_map;
    TopExp::MapShapes(shape, TopAbs_FACE, faces_map);
    Standard_Integer nFaces = faces_map.Extent();
    if (iFace < 1 || iFace > nFaces) {
        mexErrMsgIdAndTxt("StepPlotPcurves:Input","Face index out of range [1,%d]",nFaces);
    }
    const TopoDS_Face &face = TopoDS::Face(faces_map(iFace));
    TopExp_Explorer edges(face,TopAbs_EDGE);
    for (int iEdges = 0; edges.More(); edges.Next(), iEdges++) {
        const TopoDS_Edge &edge = TopoDS::Edge(edges.Current());
        plotEdgeOnFace(face, iEdges, edge);
    }
}
