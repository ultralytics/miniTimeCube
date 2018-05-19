//purpose:  photon intercepts on a rectangular surface
//inputs:   AT-attenuationTime(ns)
//outputs:  
//example:  [p2,dt,hitplane]=cubeinterceptc(p1,vel,Lr);
    
#include <math.h>
#include "mex.h"
#include <iostream>

void diffusevecc(double nx, double ny, double nz, double *X, double *Y, double *Z){
    float u,theta,b, x,y,z, yaw,pitch, cp,sp, cy,sy; 
        //u = (float)rand()/RAND_MAX; //ISOTROPIC
        u = acos((float)rand()/RAND_MAX)/(M_PI/2); //DIFFUSE, Lambert's Cosine Law https://en.wikipedia.org/wiki/Lambert%27s_cosine_law
        
        theta=(float)rand()/RAND_MAX*2*M_PI;  b=sqrt(1-u*u);
        z=b*cos(theta);  y=b*sin(theta);  x=u; //z and x switched

        yaw = atan2(ny,nx);  //pitch = asin(-nz); //unit vector! otherwise  asin(-nz/r)
        cp=sqrt(1 - nz*nz);     sp=-nz;        //cp=cos(pitch);        sp=sin(pitch);
        cy=cos(yaw);            sy=sin(yaw);
        
        *X = cp*cy*x - sy*y + cy*sp*z;
        *Y = cp*sy*x + cy*y + sp*sy*z;
        *Z = cp*z - sp*x;
    return;
}

void isovecs(double* X, double* Y, double* Z){
    float u,theta,b;
    
    u = (float)rand()/RAND_MAX*2 - 1;
    theta = (float)rand()/RAND_MAX*2*M_PI;
    b = sqrt(1-u*u);
    
    *X = b*cos(theta);
    *Y = b*sin(theta);
    *Z = u;
    return;
}

float fresnelc(double n1, double n2, double ct){
    double k1, n1ct, n1k1, n2ct, n2k1, Rs, Rp, R;
    if (ct>0){
            k1 = ( sqrt(1-pow(n1/n2,2)*(1-ct*ct)) );
            n1ct = n1*ct;  n1k1 = n1*k1;  n2ct = n2*ct;  n2k1 = n2*k1;
            Rs = pow( (n1ct - n2k1)/(n1ct + n2k1),2);
            Rp = pow( (n1k1 - n2ct)/(n1k1 + n2ct),2);
            R = fmin( (Rs+Rp)/2, 1);}
    return R;
}

void intercept1(double* P1, double* T1, double* V, int* wl, int shapeID, double* Lr, double* LT, double* RT, double* RWL,
        double* P2, double* V2, double* T2, int* S, double* WLb, double* T2b, int nX) //S=status
{
    int i, j, k, l, wli, nT=1919;
    bool sx,sy,sz;
    double t, r, tx,ty,tz, ux,uy,uz, vx,vy,vz, p1x,p1y,p1z, nx,ny,nz, speed, clight=299.792458, RN; //c=300mm/ns;
    
    for (i=0; i<nX; i++){
        j=i+nX;  k=j+nX;  l=k+nX;  wli=int(wl[i]-81);

        t = -LT[wli+2*nT]*log((float)rand()/RAND_MAX); //(ns) random attenuation time from mean attenuation time
        speed = LT[wli+5*nT]; // photon group speed (mm/ns)

        
        ux=V[i]; uy=V[j]; uz=V[k]; if (ux==0 && uy==0 && uz==0){ isovecs(&ux,&uy,&uz); }

        nx=0;  vx=ux*speed;  p1x=P1[i]-Lr[3];
        ny=0;  vy=uy*speed;  p1y=P1[j]-Lr[4];
        nz=0;  vz=uz*speed;  p1z=P1[k]-Lr[5];
        //std::cout << "r =  " << P1[i] << ",  " << Lr[3] << ",  " << p1x <<"\n";
        
        sx=vx>0;   if (sx) {tx=(Lr[0]-p1x)/vx;} else {tx=(-Lr[0]-p1x)/vx;}
        sy=vy>0;   if (sy) {ty=(Lr[1]-p1y)/vy;} else {ty=(-Lr[1]-p1y)/vy;}
        sz=vz>0;   if (sz) {tz=(Lr[2]-p1z)/vz;} else {tz=(-Lr[2]-p1z)/vz;}
        
        if (shapeID==2){ //cube
            if         (tx<ty && tx<tz && tx<t)   {t=tx; S[i]=1; if(sx){nx=-1;} else {nx=1;}}
            else if    (ty<tx && ty<tz && ty<t)   {t=ty; S[i]=2; if(sy){ny=-1;} else {ny=1;}}
            else if    (tz<t)                     {t=tz; S[i]=3; if(sz){nz=-1;} else {nz=1;}}
            P2[i]=p1x+t*vx+Lr[3];  P2[j]=p1y+t*vy+Lr[4];  P2[k]=p1z+t*vz+Lr[5];  T2[i]=T1[i]+t;
        }
        else if(shapeID==3){ //vertical cylinder
            double a,b,c,tr;
            r = Lr[0];            //std::cout << "r=" << r << "\n";
            a = vx*vx + vy*vy; //sum(velxy.^2,2);
            b = 2*(p1x*vx + p1y*vy); //2*sum(p1xy.*velxy,2);
            c = (p1x*p1x + p1y*p1y) - r*r; //sum(p1xy.^2,2) - r^2;
            tr = (-b+sqrt(b*b-4*a*c))/(2*a); //quadratic time to hit radius
            
            if         (tr<tz && tr<t)   {t=tr; S[i]=1;} //circular part
            else if    (tz<t)            {t=tz; S[i]=3; if(sz){nz=-1;} else {nz=1;}}//S0 ATTENUATED
            P2[i]=p1x+t*vx+Lr[3];  P2[j]=p1y+t*vy+Lr[4];  P2[k]=p1z+t*vz+Lr[5];  T2[i]=T1[i]+t;   if (S[i]==1){nx=-(P2[i]-Lr[3])/r; ny=-(P2[j]-Lr[4])/r;}
        }
        
        RN = (float)rand()/RAND_MAX;
        if (S[i]==0){
            if (RN < LT[wli+3*nT]){
                S[i]=4;
                T2b[i]=T2[i]+RT[(int)((float)rand()/RAND_MAX*16384)];
                WLb[i]=wl[i]+RWL[(int)((float)rand()/RAND_MAX*16384)];
            }
            continue;
        } //S0 ATTENUATED, S4 RE-EMISSION
        
        if (RN < .00001){ diffusevecc(nx,ny,nz,&V2[i],&V2[j],&V2[k]); S[i]=5; WLb[i]=wl[i];  T2b[i]=T2[i]; continue;}//5% DIFFUSE REFLECTION
        
        double ct, R;
        ct = -nx*ux - ny*uy - nz*uz;
        R = fresnelc(LT[wli], LT[wli+(5+S[i])*nT], ct); //LT5-X, LT6-Y, LT7-Z
        if (RN < R){ //SPECULAR REFLECTION
            V2[i] = ux + 2*ct*nx;
            V2[j] = uy + 2*ct*ny;
            V2[k] = uz + 2*ct*nz;  S[i]=5;  WLb[i]=wl[i];  T2b[i]=T2[i];}
        else {S[i] += 10;} //TRANSMISSION!  S11=XTRANSMISSION, S12=YTRANSMISSION, S13=ZTRANSMISSION

    }
    return;
}

// Main function definitions ----------------------------------------------
//#define IN_A		prhs[0]
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *A, *B, *C, *E, *F; //outputs
    int  div[2], nd, n0, *D; //number of dimensions and dimension vector
    const int *din;
    
    //mexPrintf("function has %d inputs and %d outputs\n",nlhs,nrhs);
    //mexPrintf("A is %dD, MxN=%dx%d\n",mxGetNumberOfDimensions(prhs[0]),mxGetM(prhs[0]),mxGetN(prhs[0]));
    //if (nrhs<5)  mexErrMsgTxt("too few inputs");
    //if (!mxIsDouble(plhs[0]) || !mxIsDouble(plhs[1])) mexErrMsgTxt("rangec: inputs must be double");
    
    //nX = fmax(mxGetM(prhs[0]),mxGetN(prhs[0]));
    n0 = mxGetM(prhs[0]); //nrows
    
    // allocate outputs, then assign pointers to outputs
    nd = mxGetNumberOfDimensions(prhs[0]);
    din = mxGetDimensions(prhs[0]);
    div[0]=din[0]; div[1]=1;
    plhs[0] = mxCreateNumericArray(nd, din, mxDOUBLE_CLASS, mxREAL);      A = mxGetPr(plhs[0]);
    plhs[1] = mxCreateNumericArray(nd, din, mxDOUBLE_CLASS, mxREAL);      B = mxGetPr(plhs[1]);
    plhs[2] = mxCreateNumericArray(nd, div, mxDOUBLE_CLASS, mxREAL);      C = mxGetPr(plhs[2]);
    plhs[3] = mxCreateNumericArray(nd, div, mxUINT32_CLASS, mxREAL);      D = (int*)mxGetPr(plhs[3]);
    plhs[4] = mxCreateNumericArray(nd, div, mxDOUBLE_CLASS, mxREAL);      E = mxGetPr(plhs[4]);
    plhs[5] = mxCreateNumericArray(nd, div, mxDOUBLE_CLASS, mxREAL);      F = mxGetPr(plhs[5]);

    
    //Do the actual computations in a subroutine
    intercept1( mxGetPr(prhs[0]), mxGetPr(prhs[1]), mxGetPr(prhs[2]), (int*) mxGetPr(prhs[3]), int(*mxGetPr(prhs[4])), mxGetPr(prhs[5]), mxGetPr(prhs[6]), mxGetPr(prhs[7]), mxGetPr(prhs[8]), 
            A, B, C, D, E, F, n0);
    
    return;
}