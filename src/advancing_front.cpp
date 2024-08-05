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
#include <functional>
#include "mex.h"
#include "matrix.h"
#include <Eigen/Dense>
#include "advancing_front.hpp"

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

template <int d, int n>
struct FunctionHandlesLoader {
    const mxArray *G_mxArray;
    const mxArray *dG_mxArray;
    const mxArray *h_mxArray;
    
    FunctionHandlesLoader(const mxArray *_G_mxArray, const mxArray *_dG_mxArray, const mxArray *_h_mxArray)
    : G_mxArray(_G_mxArray), dG_mxArray(_dG_mxArray), h_mxArray(_h_mxArray)
    {
        if (mxGetClassID(G_mxArray) != mxFUNCTION_CLASS) {
            mexErrMsgIdAndTxt("AdvancingFront:Input","G is not a function handle");
        }
        if (mxGetClassID(dG_mxArray) != mxFUNCTION_CLASS) {
            mexErrMsgIdAndTxt("AdvancingFront:Input","dG is not a function handle");
        }
        if (mxGetClassID(h_mxArray) != mxFUNCTION_CLASS) {
            mexErrMsgIdAndTxt("AdvancingFront:Input","h is not a function handle");
        }
    }
    
    Point<n> G(const Point<d> &p) const
    {
        mxArray *plhs_feval[1];
        mxArray *prhs_feval[2];
        prhs_feval[0] = (mxArray*)G_mxArray;
        prhs_feval[1] = mxArray_from_point<d>(p);
        int errcode = mexCallMATLAB(1, plhs_feval, 2, prhs_feval, "feval");
        if (errcode) {
            mexErrMsgIdAndTxt("AdvancingFront:Feval",
                              "Evaluation of G produced errors");
        }
        if (!mxIsDouble(plhs_feval[0]) || mxGetNumberOfElements(plhs_feval[0]) != n) {
            mexErrMsgIdAndTxt("AdvancingFront:Feval",
                              "Output of G must have %d elements",n);
        }
        double *Gp_buf = mxGetDoubles(plhs_feval[0]);
        Point<n> Gp;
        std::memcpy(Gp.data(), Gp_buf, n*sizeof(double));
        mxDestroyArray(plhs_feval[0]);
        mxDestroyArray(prhs_feval[1]);
        return Gp;
    }
    
    Eigen::Matrix<double,n,d> dG(const Point<d> &p) const
    {
        mxArray *plhs_feval[1];
        mxArray *prhs_feval[2];
        prhs_feval[0] = (mxArray*)dG_mxArray;
        prhs_feval[1] = mxArray_from_point<d>(p);
        int errcode = mexCallMATLAB(1, plhs_feval, 2, prhs_feval, "feval");
        if (errcode) {
            mexErrMsgIdAndTxt("AdvancingFront:Feval",
                              "Evaluation of dG produced errors");
        }
        if (!mxIsDouble(plhs_feval[0]) || mxGetNumberOfElements(plhs_feval[0]) != n*d) {
            mexErrMsgIdAndTxt("AdvancingFront:Feval",
                              "Output of dG must have %d elements",n*d);
        }
        double *dGp_buf = mxGetDoubles(plhs_feval[0]);
        Eigen::Matrix<double,n,d> dGp;
        std::memcpy(dGp.data(), dGp_buf, n*d*sizeof(double));
        mxDestroyArray(plhs_feval[0]);
        mxDestroyArray(prhs_feval[1]);
        return dGp;
    }
    
    double h(const Point<n> &x) const
    {
        mxArray *plhs_feval[1];
        mxArray *prhs_feval[2];
        prhs_feval[0] = (mxArray*)h_mxArray;
        prhs_feval[1] = mxArray_from_point<n>(x);
        int errcode = mexCallMATLAB(1, plhs_feval, 2, prhs_feval, "feval");
        if (errcode) {
            mexErrMsgIdAndTxt("AdvancingFront:Feval",
                              "Evaluation of h produced errors");
        }
        if (!mxIsDouble(plhs_feval[0]) || mxGetNumberOfElements(plhs_feval[0]) != 1) {
            mexErrMsgIdAndTxt("AdvancingFront:Feval",
                              "Output of h must be a scalar");
        }
        double hx = mxGetScalar(plhs_feval[0]);
        mxDestroyArray(plhs_feval[0]);
        mxDestroyArray(prhs_feval[1]);
        return hx;
    }
};

template <int d, int n>
void mexFunctionTemplated(mxArray *plhs[], const mxArray *prhs[])
// function [Y,GY,nGY] = advancing_front(aabb,Z_points,Z_normals,Y0,G,dG,h,Nmax,seed)
{
    // Rename input and output arguments
    const mxArray *aabb_mxArray = prhs[0];
    const mxArray *Z_points_mxArray = prhs[1];
    const mxArray *Z_normals_mxArray = prhs[2];
    const mxArray *Y0_mxArray = prhs[3];
    const mxArray *G_mxArray = prhs[4];
    const mxArray *dG_mxArray = prhs[5];
    const mxArray *h_mxArray = prhs[6];
    const mxArray *Nmax_mxArray = prhs[7];
    const mxArray *seed_mxArray = prhs[8];
    mxArray *Y_mxArray = nullptr;
    mxArray *GY_mxArray = nullptr;
    mxArray *nGY_mxArray = nullptr;
    
    // Unwrap aabb
    double *aabb_buf = mxGetDoubles(aabb_mxArray);
    if (aabb_buf == nullptr) {
        mexErrMsgIdAndTxt("AdvancingFront:Input","aabb has the wrong type");
    }
    if (mxGetM(aabb_mxArray) != 2 || mxGetN(aabb_mxArray) != (size_t)d) {
        mexErrMsgIdAndTxt("AdvancingFront:Input","aabb must have %d rows and %d columns",2,d);
    }
    Point<d> aabb_min;
    Point<d> aabb_max;
    for (size_t i = 0; i < d; i++) {
        aabb_min(i) = aabb_buf[2*i];
        aabb_max(i) = aabb_buf[2*i+1];
    }
    AABB<d> aabb(aabb_min, aabb_max);
    
    // Unwrap Z_points and Z_normals
    if (!mxIsEmpty(Z_points_mxArray) && !mxIsDouble(Z_points_mxArray)) {
        mexErrMsgIdAndTxt("AdvancingFront:Input","Z_points must have type double");
    }
    if (!mxIsEmpty(Z_normals_mxArray) && !mxIsDouble(Z_normals_mxArray)) {
        mexErrMsgIdAndTxt("AdvancingFront:Input","Z_normals must have type double");
    }
    size_t M_Z_points = mxGetM(Z_points_mxArray);
    size_t N_Z_points = mxGetN(Z_points_mxArray);
    size_t M_Z_normals = mxGetM(Z_normals_mxArray);
    size_t N_Z_normals = mxGetN(Z_normals_mxArray);
    if (M_Z_points != M_Z_normals) {
        mexErrMsgIdAndTxt("AdvancingFront:Input",
                          "Z_points and Z_normals must have the same number of rows");
    }
    if (N_Z_points != N_Z_normals) {
        mexErrMsgIdAndTxt("AdvancingFront:Input",
                          "Z_points and Z_normals must have the same number of columns");
    }
    if (N_Z_points != 0 && N_Z_points != d) {
        mexErrMsgIdAndTxt("AdvancingFront:Input",
                          "Z_points must have either 0 or %d columns",d);
    }
    double *Z_points_buf = mxGetDoubles(Z_points_mxArray);
    double *Z_normals_buf = mxGetDoubles(Z_normals_mxArray);
    Nodes<d> Z(M_Z_points, Z_points_buf, Z_normals_buf);
    
    // Unwrap Y0
    if (!mxIsEmpty(Y0_mxArray) && !mxIsDouble(Y0_mxArray)) {
        mexErrMsgIdAndTxt("AdvancingFront:Input","Y0 must have type double");
    }
    size_t M_Y0 = mxGetM(Y0_mxArray);
    size_t N_Y0 = mxGetN(Y0_mxArray);
    if (N_Y0 != 0 && N_Y0 != d) {
        mexErrMsgIdAndTxt("AdvancingFront:Input","Y0 must have either 0 or %d columns",d);
    }
    double *Y0_buf = mxGetDoubles(Y0_mxArray);
    Nodes<d> Y(M_Y0, Y0_buf);
    
    // Unwrap G, dG, h
    FunctionHandlesLoader<d,n> fhl(G_mxArray,dG_mxArray,h_mxArray);
    
    // Unwrap Nmax
    if (!mxIsScalar(Nmax_mxArray) || !mxIsNumeric(Nmax_mxArray)) {
        mexErrMsgIdAndTxt("AdvancingFront:Input","Nmax must be a scalar of numeric type");
    }
    size_t Nmax = (size_t)mxGetScalar(Nmax_mxArray);

    // Unwrap seed
    if (!mxIsScalar(seed_mxArray) || !mxIsNumeric(seed_mxArray)) {
        mexErrMsgIdAndTxt("AdvancingFront:Input","seed must be a scalar of numeric type");
    }
    uint64_t seed = (uint64_t)mxGetScalar(seed_mxArray);
    std::mt19937_64 rng(seed);
    
    Nodes<n> GY = advancing_front<d,n>(aabb, Z, Y,
        std::bind(&FunctionHandlesLoader<d,n>::G,fhl,std::placeholders::_1),
        std::bind(&FunctionHandlesLoader<d,n>::dG,fhl,std::placeholders::_1),
        std::bind(&FunctionHandlesLoader<d,n>::h,fhl,std::placeholders::_1), Nmax, rng);
    
    // Create Y
    size_t NY = Y.points.size();
    Y_mxArray = mxCreateDoubleMatrix(NY, d, mxREAL);
    double *Y_buf = mxGetDoubles(Y_mxArray);
    for (size_t i = 0; i < NY; i++) {
        for (size_t j = 0; j < d; j++) {
            Y_buf[i+j*NY] = Y.points[i][j];
        }
    }
    
    // Create GY
    size_t NGY = GY.points.size();
    GY_mxArray = mxCreateDoubleMatrix(NGY, n, mxREAL);
    double *GY_buf = mxGetDoubles(GY_mxArray);
    for (size_t i = 0; i < NGY; i++) {
        for (size_t j = 0; j < n; j++) {
            GY_buf[i+j*NY] = GY.points[i][j];
        }
    }
    
    // Create nGY
    size_t NnGY = GY.normals.size();
    nGY_mxArray = mxCreateDoubleMatrix(NnGY, n, mxREAL);
    double *nGY_buf = mxGetDoubles(nGY_mxArray);
    for (size_t i = 0; i < NnGY; i++) {
        for (size_t j = 0; j < n; j++) {
            nGY_buf[i+j*NY] = GY.normals[i][j];
        }
    }
    
    plhs[0] = Y_mxArray;
    plhs[1] = GY_mxArray;
    plhs[2] = nGY_mxArray;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
// function [Y,GY,nGY] = advancing_front(d,n,aabb,Z_points,Z_normals,Y0,G,dG,h,Nmax,seed)
{
    // Rename input and output arguments
    if (nrhs != 11) {
        mexErrMsgIdAndTxt("AdvancingFront:Input","Wrong number of input arguments");
    }
    if (nlhs != 3) {
        mexErrMsgIdAndTxt("AdvancingFront:Input","Wrong number of output arguments");
    }
    const mxArray *d_mxArray = prhs[0];
    const mxArray *n_mxArray = prhs[1];
    
    // Unwrap d
    if (!mxIsScalar(d_mxArray) || !mxIsNumeric(d_mxArray)) {
        mexErrMsgIdAndTxt("AdvancingFront:Input","d must be a scalar of numeric type");
    }
    int d = (int)mxGetScalar(d_mxArray);
    if (d < 1 || d > 3) {
        mexErrMsgIdAndTxt("AdvancingFront:Input","d must be in the range [1,3]");
    }
    
    // Unwrap n
    if (!mxIsScalar(n_mxArray) || !mxIsNumeric(n_mxArray)) {
        mexErrMsgIdAndTxt("AdvancingFront:Input","n must be a scalar of numeric type");
    }
    int n = (int)mxGetScalar(n_mxArray);
    if (n < d || n > 3) {
        mexErrMsgIdAndTxt("AdvancingFront:Input","n must be in the range [d,3]");
    }
    
    // Call the appropriate templated mexFunction
    if (d == 1 && n == 1) {
        mexFunctionTemplated<1,1>(plhs,&prhs[2]);
    } else if (d == 1 && n == 2) {
        mexFunctionTemplated<1,2>(plhs,&prhs[2]);
    } else if (d == 1 && n == 3) {
        mexFunctionTemplated<1,3>(plhs,&prhs[2]);
    } else if (d == 2 && n == 2) {
        mexFunctionTemplated<2,2>(plhs,&prhs[2]);
    } else if (d == 2 && n == 3) {
        mexFunctionTemplated<2,3>(plhs,&prhs[2]);
    } else if (d == 3 && n == 3) {
        mexFunctionTemplated<3,3>(plhs,&prhs[2]);
    } else {
        mexErrMsgIdAndTxt("AdvancingFront:Input","Invalid pair (d,n)");
    }
}
