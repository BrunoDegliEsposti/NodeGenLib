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
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Shape.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
// function [] = step_print_faces_info(filename)
{
    // Rename input arguments
    if (nrhs != 1) {
        mexErrMsgIdAndTxt("StepPrintFacesInfo:Input","Function takes in only one argument");
    }
    const mxArray *filename_mxArray = prhs[0];
    
    // Unwrap filename
    Standard_CString filename = mxArrayToString(filename_mxArray);
    if (filename == nullptr) {
        mexErrMsgIdAndTxt("StepPrintFacesInfo:Input","Filename must be a char array");
    }

    // Read input file
    STEPControl_Reader reader;
    IFSelect_ReturnStatus read_status = reader.ReadFile(filename);
    if (read_status != IFSelect_RetDone) {
        mexErrMsgIdAndTxt("StepPrintFacesInfo:Input","Error reading STEP file");
    }
    Standard_Integer nShapes = reader.TransferRoots();
    TopoDS_Shape shape = reader.OneShape();

    // Iterate over all unique faces
    TopTools_IndexedMapOfShape faces_map;
    TopExp::MapShapes(shape, TopAbs_FACE, faces_map);
    Standard_Integer nFaces = faces_map.Extent();
    for (int iFace = 1; iFace <= nFaces; iFace++) {
        const TopoDS_Face &face = TopoDS::Face(faces_map(iFace));
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
        mexPrintf("  Bounding box of pcurves: (%f,%f) x (%f,%f)\n",u0,u1,v0,v1);
    }
}
