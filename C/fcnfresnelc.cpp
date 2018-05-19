//purpose:  Finds the bsx ranges between two sets of 2D or 3D vectors A and B
//inputs:   A (nx3), optional B (mx3)
//outputs:  range (nxm), range^2 (nxm), dx (nxm), dy (nxm), dz (nxm)
//example:  R = fcnfresnelc(n1,n2,ct)

#include <math.h>
#include "mex.h"
void cfunction(double* n_1, double* n_2, double* ct, double* R, int MN1, int MN3)
{
    int i;
    double n1, n2, a, k1, n1ct, n1k1, n2ct, n2k1, Rs, Rp;
    
    
    if (MN1==1){ //only one set of n1 and n2 supplied, but maybe many ct
        n1 = n_1[0];
        n2 = n_2[0];
        a = pow(n1/n2,2);
        for (i=0; i<MN3; i++){
            if (ct[i]>0){
                k1 = ( sqrt(1-a*(1-pow(ct[i],2))) );
                n1ct = n1*ct[i];
                n1k1 = n1*k1;
                n2ct = n2*ct[i];
                n2k1 = n2*k1;
                Rs = pow( (n1ct - n2k1)/(n1ct + n2k1),2);
                Rp = pow( (n1k1 - n2ct)/(n1k1 + n2ct),2);
                R[i] = fmin( (Rs+Rp)/2, 1);
            }
        }}
    else{
        for (i=0; i<MN3; i++){
            if (ct[i]>0){
                n1 = n_1[i];
                n2 = n_2[i];
                a = pow(n1/n2,2);
                k1 = ( sqrt(1-a*(1-pow(ct[i],2))) );
                n1ct = n1*ct[i];
                n1k1 = n1*k1;
                n2ct = n2*ct[i];
                n2k1 = n2*k1;
                Rs = pow( (n1ct - n2k1)/(n1ct + n2k1),2);
                Rp = pow( (n1k1 - n2ct)/(n1k1 + n2ct),2);
                R[i] = fmin( (Rs+Rp)/2, 1);
            }}
    }
    return;
}

// Main function definitions ----------------------------------------------
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *A; //outputs
    int  MN1, MN3, nd, nX; //number of dimensions and dimension vector
    const int *din;
    
    //mexPrintf("function has %d inputs and %d outputs\n",nlhs,nrhs);
    //mexPrintf("A is %dD, MxN=%dx%d\n",mxGetNumberOfDimensions(prhs[0]),mxGetM(prhs[0]),mxGetN(prhs[0]));
    if (nlhs>1)  mexErrMsgTxt("too many outputs");
    if (nrhs<3)  mexErrMsgTxt("too few inputs");
    //if (!mxIsDouble(plhs[0]) || !mxIsDouble(plhs[1])) mexErrMsgTxt("rangec: inputs must be double");
    
    MN1 = mxGetM(prhs[0])*mxGetN(prhs[0]);
    MN3 = mxGetM(prhs[2])*mxGetN(prhs[2]);
    
    // allocate outputs, then assign pointers to outputs
    nd = mxGetNumberOfDimensions(prhs[2]);
    din = mxGetDimensions(prhs[2]);    
    plhs[0] = mxCreateNumericArray(nd, din, mxDOUBLE_CLASS, mxREAL);      A = mxGetPr(plhs[0]);

    // Do the actual computations in a subroutine
    cfunction(mxGetPr(prhs[0]),mxGetPr(prhs[1]),mxGetPr(prhs[2]), A, MN1, MN3);
    
    return;
}