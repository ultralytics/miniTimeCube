//purpose:  Finds the bsx ranges between two sets of 2D or 3D vectors A and B
//inputs:   A (nx3), optional B (mx3)
//outputs:  range (nxm), range^2 (nxm), dx (nxm), dy (nxm), dz (nxm)
//example:  [r,rs,dx,dy,dz] = rangec(A,B)

#include <vector>
#include <iostream>
#include <math.h>
#include <cmath>
#include "mex.h"

mxArray * getMexArrayN1(const std::vector<double>& v){mxArray * mx = mxCreateDoubleMatrix(v.size(), 1, mxREAL);  std::copy(v.begin(), v.end(), mxGetPr(mx));  return mx;}
mxArray * getMexArrayN3(const std::vector<double>& v){mxArray * mx = mxCreateDoubleMatrix(v.size()/3,  3, mxREAL);  std::copy(v.begin(), v.end(), mxGetPr(mx));  return mx;}
mxArray * getMexArrayNM(const std::vector<double>& v, int rows, int cols){mxArray * mx = mxCreateDoubleMatrix(rows,  cols, mxREAL);  std::copy(v.begin(), v.end(), mxGetPr(mx));  return mx;}
std::vector<double> vecInput(const mxArray *a){double *A=mxGetPr(a);  std::vector<double> w(A, A + mxGetM(a)*mxGetN(a) ); return w;}


void multiplyBy2(std::vector<double>& v){
     v[1]=16;
    return;
}


void mexFunction(int nlhs, mxArray *plhs[ ], int nrhs, const mxArray *prhs[ ]) {
    //double *A=mxGetPr(prhs[0]);      
    //std::vector<double> w(A, A + mxGetM(prhs[0])*mxGetN(prhs[0]) ); // copy A to w
    //std::vector<double> v;  v.resize(12);  std::copy(A, A + 5, v.begin()); // copy A to v

//     v.push_back(0);
    
    std::vector<double> w = vecInput(prhs[0]);
    
    multiplyBy2(w);

    plhs[0] = getMexArrayN3(w);
}