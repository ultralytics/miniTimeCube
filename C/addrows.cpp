//purpose:  Finds the bsx ranges between two sets of 2D or 3D vectors A and B
//inputs:   A (nx3), optional B (mx3)
//outputs:  range (nxm), range^2 (nxm), dx (nxm), dy (nxm), dz (nxm)
//example:  [r,rs,dx,dy,dz] = rangec(A,B)

#include <vector>
#include <iostream>
#include <math.h>
#include <cmath>
#include "mex.h"



std::vector<double> example(std::vector<double> Xa){
    std::vector<double> Xb;
    Xb.resize(20);
    for(int i=0; i<Xa.size(); i++)
    {
        Xb[i]=i*2;
    }
    Xb.resize(10); //trim excess rows
    return Xb;
}


void cfunction(int rows, double* A)
{
    int i;
    float u,theta,b;
    
    for (i=0; i<(rows+1); i++)
    {
        u = (float)rand()/RAND_MAX*2 - 1;
        theta = (float)rand()/RAND_MAX*2*M_PI;
        b = sqrt(1-u*u);
        
        A[i] = b*cos(theta);
        A[i+rows] = b*sin(theta);
        A[i+2*rows] = u;
    }
     return;
}


void cfunction2(int rows, std::vector<double> A)
{
    int i;
    float u,theta,b;
    
    for (i=0; i<(rows+1); i++)
    {
        u = (float)rand()/RAND_MAX*2 - 1;
        theta = (float)rand()/RAND_MAX*2*M_PI;
        b = sqrt(1-u*u);
        
        A[i] = b*cos(theta);
        A[i+rows] = b*sin(theta);
        A[i+2*rows] = u;
    }
     return;
}



// Main function definitions ----------------------------------------------
// #define IN_A		prhs[0]
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *A, *B; //outputs
    int  nd, *d, nrows; //number of dimensions and dimension vector

    //mexPrintf("function has %d inputs and %d outputs\n",nlhs,nrhs);
    //mexPrintf("A is %dD, MxN=%dx%d\n",mxGetNumberOfDimensions(prhs[0]),mxGetM(prhs[0]),mxGetN(prhs[0]));
    if (nlhs>1)  mexErrMsgTxt("too many outputs");
    if (nrhs<1)  mexErrMsgTxt("too few inputs");
    //if (!mxIsDouble(plhs[0]) || !mxIsDouble(plhs[1])) mexErrMsgTxt("rangec: inputs must be double");
    
    nrows = mxGetScalar(prhs[0]);
    //mexPrintf("function has %d rows\n",nrows);

    // allocate outputs, then assign pointers to outputs
    nd = 2;  d = (int*) calloc(nd, sizeof(int));  d[0]=nrows;  d[1]=3;    
    plhs[0] = mxCreateNumericArray(nd, d, mxDOUBLE_CLASS, mxREAL);      A = mxGetPr(plhs[0]);

    
    std::vector<double> Xb;
    Xb.resize(16);
    cfunction2(nrows, Xb);
    std::cout << "r =  " << Xb[0] << ",  " << Xb[1] << ",  " << Xb[2] <<"\n";

    
    // Do the actual computations in a subroutine
    //cfunction(nrows, A);
    
    
    
    
    
    //cfunction3(nrows);
    
    return;
}