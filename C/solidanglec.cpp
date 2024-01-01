//purpose:  Finds the bsx ranges between two sets of 2D or 3D vectors A and B
//inputs:   A (nx3), optional B (mx3)
//outputs:  range (nxm), range^2 (nxm), dx (nxm), dy (nxm), dz (nxm)
//example:  [r,rs,dx,dy,dz] = rangec(A,B)

#include <vector>
#include <iostream>
#include <math.h>
#include <cmath>
#include "mex.h"

float interpFt(double* X, double XI){
    int i = (int) (XI * 2999); //X is Ft lookup table, 1x3000 from ct=0 to 1
    //std::cout << "ct=" << XI << ", index is " << i << "\n";
    return(X[i]);
}

std::vector<double> fcngetreflections(std::vector<double> Xa){
    std::vector<double> Xb;
    Xb.resize(20);
    for(int i=0; i<Xa.size(); i++)
    {
        
    }
    Xb.resize(10); //trim excess rows
    return Xb;
}


void cfunction(double* X, double* P, double* N, double* pa, double* al, double* kFty, double* Fall, double* r, int MX, int MP)
//X=pointxyz, P=pixelxyz, N=pixelnormalvecs, pa=pixelarea(mm^2), al=attenuation length (mm)
{
    int i, j, k;
    double Px,Py,Pz,dx,dy,dz,cx,cy,cz,x,y,z,rs,dotproduct, Fsa, Fna, Ft, fr, ct, ctx, nial=-1/al[0], res=pa[0]/M_PI; //nial=negative inverse attenuation length, res = pixel radius equivalent squared
    std::vector<double> Xa;  Xa.resize(7); //Xa = [x y z nxr nyr nzr gen]

    for (i=0; i<MX; i++) //X (point sources), nx6, [xyz,mx,my,mz], mx=number of reflectiosn through x plane etc.
    {
        cx = X[i+3*MX];  x=X[i];
        cy = X[i+4*MX];  y=X[i+MX];
        cz = X[i+5*MX];  z=X[i+2*MX];
        
        for (j=0; j<MP; j++) //P (pixels)
        {
            k = j+i*MP;
            dx = P[j] - x;
            dy = P[j+MP] - y;
            dz = P[j+2*MP] - z;
            dotproduct = N[j]*dx + N[j+MP]*dy + N[j+2*MP]*dz;
            
            if (dotproduct>0) // possible transmission
            {            
                rs = dx*dx + dy*dy + dz*dz;
                r[k] = sqrt(rs);
                ct = dotproduct/r[k]; //cos(theta) = N, (X-P) dot product divided by range
                Ft = interpFt(kFty, ct);
                if (Ft>0)
                {
                    if (cx>0) {ctx=std::abs(dx)/r[k];  fr=1-interpFt(kFty,ctx);  if(cx>1){fr=pow(fr,cx);};  Ft*=fr;}
                    if (cy>0) {ctx=std::abs(dy)/r[k];  fr=1-interpFt(kFty,ctx);  if(cy>1){fr=pow(fr,cy);};  Ft*=fr;}
                    if (cz>0) {ctx=std::abs(dz)/r[k];  fr=1-interpFt(kFty,ctx);  if(cz>1){fr=pow(fr,cz);};  Ft*=fr;}

                    Fsa = 0.5*(1 - r[k]/sqrt(rs+res*ct));
                    Fna = exp(r[k]*nial);
                    Fall[k] = Fsa * Fna * Ft;
                }
            } //if (dotproduct>0)
        }//for (i=0; i<MX; i++)
    }
    return;
}


// Main function definitions ----------------------------------------------
//#define IN_A		prhs[0]
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *A, *B; //outputs
    int  nd, *d; //number of dimensions and dimension vector
    
    //mexPrintf("function has %d inputs and %d outputs\n",nlhs,nrhs);
    //mexPrintf("A is %dD, MxN=%dx%d\n",mxGetNumberOfDimensions(prhs[0]),mxGetM(prhs[0]),mxGetN(prhs[0]));
    if (nlhs>6)  mexErrMsgTxt("too many outputs");
    if (nrhs<2)  mexErrMsgTxt("too few inputs");
    //if (!mxIsDouble(plhs[0]) || !mxIsDouble(plhs[1])) mexErrMsgTxt("inputs must be double");

    // allocate outputs, then assign pointers to outputs
    nd = 2;  d = (int*) calloc(nd, sizeof(int));  *d=mxGetM(prhs[1]);  d[1]=mxGetM(prhs[0]);
    plhs[0] = mxCreateNumericArray(nd, d, mxDOUBLE_CLASS, mxREAL);      A = mxGetPr(plhs[0]); //Fall
    plhs[1] = mxCreateNumericArray(nd, d, mxDOUBLE_CLASS, mxREAL);      B = mxGetPr(plhs[1]); //r

    // Do the actual computations in a subroutine
    cfunction(mxGetPr(prhs[0]), mxGetPr(prhs[1]), mxGetPr(prhs[2]), mxGetPr(prhs[3]), 
            mxGetPr(prhs[4]), mxGetPr(prhs[5]), A, B, mxGetM(prhs[0]), mxGetM(prhs[1]));
    
    free((void*)d);
    return;
}