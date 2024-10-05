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
#include <TopExp_Explorer.hxx>
#include <ShapeAnalysis_ShapeContents.hxx>

const char* shapeTypeToString(const TopoDS_Shape &shape)
{
    switch (shape.ShapeType()) {
        case TopAbs_COMPOUND: return "Compound";
        case TopAbs_COMPSOLID: return "Compsolid";
        case TopAbs_SOLID: return "Solid";
        case TopAbs_SHELL: return "Shell";
        case TopAbs_FACE: return "Face";
        case TopAbs_WIRE: return "Wire";
        case TopAbs_EDGE: return "Edge";
        case TopAbs_VERTEX: return "Vertex";
        default:
            mexErrMsgIdAndTxt("StepPrintContents:Input","Unknown shape type");
            return nullptr;
    }
}

void printShapeTypeTree(const TopoDS_Shape &shape, int indentationLevel)
{
    for (int i = 0; i < indentationLevel*2; i++) {
        mexPrintf(" ");
    }
    mexPrintf("%s\n",shapeTypeToString(shape));
    TopAbs_ShapeEnum st = shape.ShapeType();
    if (st == TopAbs_COMPOUND || st == TopAbs_COMPSOLID || st == TopAbs_SOLID) {
        TopoDS_Iterator children(shape);
        for (; children.More(); children.Next()) {
            const TopoDS_Shape &child = children.Value();
            printShapeTypeTree(child,indentationLevel+1);
        }
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
// function [] = step_print_contents(filename)
{
    // Rename input arguments
    if (nrhs != 1) {
        mexErrMsgIdAndTxt("StepPrintContents:Input","Function takes in only one argument");
    }
    const mxArray *filename_mxArray = prhs[0];
    
    // Unwrap filename
    Standard_CString filename = mxArrayToString(filename_mxArray);
    if (filename == nullptr) {
        mexErrMsgIdAndTxt("StepPrintContents:Input","Filename must be a char array");
    }

    // Read input file
    STEPControl_Reader reader;
    IFSelect_ReturnStatus read_status = reader.ReadFile(filename);
    if (read_status != IFSelect_RetDone) {
        mexErrMsgIdAndTxt("StepPrintContents:Input","Error reading STEP file");
    }
    Standard_Integer nShapes = reader.TransferRoots();
    TopoDS_Shape shape = reader.OneShape();
    mexPrintf("Number of roots successfully translated: %d\n",nShapes);

    // Process shape
    mexPrintf("Topological hierarchy up to shells:\n");
    printShapeTypeTree(shape,0);

    // Explore shape
    ShapeAnalysis_ShapeContents contents;
    contents.Perform(shape);
    mexPrintf("Number of solids: %d unique, %d non-unique\n",
              contents.NbSharedSolids(),contents.NbSolids());
    mexPrintf("Number of shells: %d unique, %d non-unique\n",
              contents.NbSharedShells(),contents.NbShells());
    mexPrintf("Number of faces: %d unique, %d non-unique, %d free from shells\n",
              contents.NbSharedFaces(),contents.NbFaces(),contents.NbFreeFaces());
    mexPrintf("Number of wires: %d unique, %d non-unique, %d free from faces\n",
              contents.NbSharedWires(),contents.NbWires(),contents.NbFreeWires());
    mexPrintf("Number of edges: %d unique, %d non-unique, %d free from wires\n",
              contents.NbSharedEdges(),contents.NbEdges(),contents.NbFreeEdges());
    mexPrintf("Number of vertices: %d unique, %d non-unique\n",
              contents.NbSharedVertices(),contents.NbVertices());
}
