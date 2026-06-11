//purpose:  photon intercepts on a rectangular surface
//inputs:   AT-attenuationTime(ns)
//outputs:  
//example:  [p2,dt,hitplane]=cubeinterceptc(p1,vel,Lr);

#include <vector>
#include <math.h>
#include "mex.h"
#include <iostream>


mxArray * getMexArrayN1(const std::vector<double>& v){mxArray * mx = mxCreateDoubleMatrix(v.size(), 1, mxREAL);  std::copy(v.begin(), v.end(), mxGetPr(mx));  return mx;}

//mxArray * getMexArrayN1int(const std::vector<int>& v){mxArray * mx = mxCreateNumericArray(1 , v.size(), mxUINT32_CLASS, mxREAL);  std::copy(v.begin(), v.end(), mxGetPr(mx));  return mx;}

mxArray * getMexArrayN3(const std::vector<double>& v){mxArray * mx = mxCreateDoubleMatrix(v.size()/3,  3, mxREAL);  std::copy(v.begin(), v.end(), mxGetPr(mx));  return mx;}
mxArray * getMexArray3N(const std::vector<double>& v){mxArray * mx = mxCreateDoubleMatrix(3,  v.size()/3, mxREAL);  std::copy(v.begin(), v.end(), mxGetPr(mx));  return mx;}
mxArray * getMexArrayNM(const std::vector<double>& v, int rows, int cols){mxArray * mx = mxCreateDoubleMatrix(rows,  cols, mxREAL);  std::copy(v.begin(), v.end(), mxGetPr(mx));  return mx;}
std::vector<double> vecInput(const mxArray *a){double *A=mxGetPr(a);  std::vector<double> w(A, A + mxGetM(a)*mxGetN(a) ); return w;}

// void randcdfc(double* cdf, int elements, double* A){
//     *A = cdf[ (int) ((float)rand()/RAND_MAX*elements) ];
//     return;
// }

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

void fcnPhotons(std::vector<double>& P1,        //1 P1
        std::vector<double>& T1,        //2 T1
        std::vector<double>& V1,         //3 V1
        int shapeID,                    //5 ShapeID
        double* Lr,        //6 Lr
        double* LT,        //7 T
        double* RT,        //8 random Time re-emission, nR=16384
        double* RWL,       //9 random Wavelength re-emission, nR=16384
        std::vector<double>& P2,     //10 p2
        std::vector<double>& T2,     //11 t2
        std::vector<int>& S,      //12 status
        std::vector<int>& gen,    //13 generation
        std::vector<double>& wl,     //14 wavelength
        std::vector<int>& ancestor) //15 ancestor index
{
    int i, j, k, l, ia, wli, istart, iend, generation, newPhotons=T1.size(), nR=16384, nT=1919;
    bool sx,sy,sz;
    double t, r, tx,ty,tz, ux,uy,uz, vx,vy,vz, v2x,v2y,v2z, p1x,p1y,p1z, nx,ny,nz, speed, clight=299.792458, RN; //c=300mm/ns;    
    std::vector<double> parent;  parent.resize(newPhotons); //Xa = [x y z nxr nyr nzr gen]
    
     istart=0; iend=newPhotons;
     for (generation=0; generation<5000; generation++){
         //std::cout << "generation =  " << generation << "\n";
         if(newPhotons==0 || generation>5000){break;}

         newPhotons=0;
         for (i=istart; i<iend; i++){
             
             if (generation==0){
                 ancestor[i]=i;
                 parent[i]=i;
             }
            else { //inherit parent properties                
                ia = parent[i];  ancestor.push_back(ancestor[ia]);
                S.push_back(0);  T2.push_back(0);  P2.push_back(0);  P2.push_back(0);  P2.push_back(0);  gen.push_back(generation);
                P1.push_back(P2[ia*3]);  P1.push_back(P2[ia*3+1]);  P1.push_back(P2[ia*3+2]);
                if (S[ia]==4){
                    isovecs(&v2x,&v2y,&v2z); V1.push_back(v2x); V1.push_back(v2y); V1.push_back(v2z);
                    T1.push_back(T2[ia] + RT[  (int) ((float)rand()/RAND_MAX*16384) ] );
                    wl.push_back(wl[ia] + RWL[ (int) ((float)rand()/RAND_MAX*16384) ] );
                }
                else{
                    T1.push_back(T2[ia]);
                    wl.push_back(wl[ia]);
                }
            }
             
            j=i*3;  k=j+1;  l=k+1;  wli=int(wl[i]-81);
            
            t = -LT[wli+2*nT]*log((float)rand()/RAND_MAX); //(ns) random attenuation time from mean attenuation time
            speed = LT[wli+5*nT]; // photon group speed (mm/ns)
            
            ux=V1[j]; uy=V1[k]; uz=V1[l]; //if (ux==0 && uy==0 && uz==0){ isovecs(&ux,&uy,&uz); }
            
            nx=0;  vx=ux*speed;  p1x=P1[j]-Lr[3];
            ny=0;  vy=uy*speed;  p1y=P1[k]-Lr[4];
            nz=0;  vz=uz*speed;  p1z=P1[l]-Lr[5];
            
            sx=vx>0;   if (sx) {tx=(Lr[0]-p1x)/vx;} else {tx=(-Lr[0]-p1x)/vx;}
            sy=vy>0;   if (sy) {ty=(Lr[1]-p1y)/vy;} else {ty=(-Lr[1]-p1y)/vy;}
            sz=vz>0;   if (sz) {tz=(Lr[2]-p1z)/vz;} else {tz=(-Lr[2]-p1z)/vz;}
            //std::cout << "r =  " << t << "," << P1[j] << ",  " << P1[k] << ",  " << P1[l] <<"\n";


            if (shapeID==2){ //cube
                if         (tx<ty && tx<tz && tx<t)   {t=tx; S[i]=1; if(sx){nx=-1;} else {nx=1;}}
                else if    (ty<tx && ty<tz && ty<t)   {t=ty; S[i]=2; if(sy){ny=-1;} else {ny=1;}}
                else if    (tz<t)                     {t=tz; S[i]=3; if(sz){nz=-1;} else {nz=1;}}
                P2[j]=p1x+t*vx+Lr[3];  P2[k]=p1y+t*vy+Lr[4];  P2[l]=p1z+t*vz+Lr[5];  T2[i]=T1[i]+t;
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
                P2[j]=p1x+t*vx+Lr[3];  P2[k]=p1y+t*vy+Lr[4];  P2[l]=p1z+t*vz+Lr[5];  T2[i]=T1[i]+t;   if (S[i]==1){nx=-(P2[j]-Lr[3])/r; ny=-(P2[k]-Lr[4])/r;}
            }

            
            RN = (float)rand()/RAND_MAX;
            if (S[i]==0){
                if (RN < LT[wli+3*nT]){
                    S[i]=4;  parent.push_back(i);  newPhotons++;}
                continue;
            } //S0 ATTENUATED, S4 RE-EMISSION
            
            
            if (RN < .00001){
                diffusevecc(nx,ny,nz,&v2x,&v2y,&v2z);  V1.push_back(v2x); V1.push_back(v2y); V1.push_back(v2z); 
                S[i]=5;  parent.push_back(i);  newPhotons++;
                continue;
            }//5% DIFFUSE REFLECTION
            
            
            double ct, R;
            ct = -nx*ux - ny*uy - nz*uz;
            R = fresnelc(LT[wli], LT[wli+(5+S[i])*nT], ct); //LT5-X, LT6-Y, LT7-Z
            if (RN < R){ //SPECULAR REFLECTION
                v2x = ux + 2*ct*nx;
                v2y = uy + 2*ct*ny;
                v2z = uz + 2*ct*nz;  V1.push_back(v2x); V1.push_back(v2y); V1.push_back(v2z);
                S[i]=5;  parent.push_back(i);  newPhotons++;}
            else {S[i] += 10;
            } //TRANSMISSION!  S11=XTRANSMISSION, S12=YTRANSMISSION, S13=ZTRANSMISSION
             
             
         } // for photon
         istart=iend; iend=istart+newPhotons;
         
     } // for generation
    return;
}

// Main function definitions ----------------------------------------------
//#define IN_A		prhs[0]
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    std::vector<double> A=vecInput(prhs[0]);
    std::vector<double> B=vecInput(prhs[1]);
    std::vector<double> V1=vecInput(prhs[2]);
    std::vector<double> G=vecInput(prhs[3]);
    
    std::vector<double> C;     C.resize(B.size()*3);
    std::vector<double> D;     D.resize(B.size());
    std::vector<int> E;        E.resize(B.size());
    std::vector<int> F;        F.resize(B.size());
    std::vector<int> H;        H.resize(B.size());
    
    fcnPhotons(A,                      //1 P1
            B,                         //2 T1
            V1,                        //3 V1
            int(*mxGetPr(prhs[4])),  //5 ShapeID
            mxGetPr(prhs[5]),        //6 Lr
            mxGetPr(prhs[6]),        //7 T
            mxGetPr(prhs[7]),        //8 RT
            mxGetPr(prhs[8]),        //9 RWL
            C ,                         //10 P2
            D,                          //11 T2
            E,                          //12 status
            F,                          //13 generation
            G,                          //14 wavelength
            H);                         //15 parent Photon index
    
//[p1,t1,p2,t2,status,generation,wl,ia]=photontransportcV3(p1,t1,veluvecs,uint32(wl-81),shapeID,Lr,T,RT,RWL);
    std::vector<double> EV(E.begin(), E.end());
    std::vector<double> FV(F.begin(), F.end());
    std::vector<double> HV(H.begin(), H.end());
    
    plhs[0] = getMexArray3N(A); //P1
    plhs[1] = getMexArrayN1(B); //T1
    plhs[2] = getMexArray3N(C); //P2
    plhs[3] = getMexArrayN1(D); //T2
    plhs[4] = getMexArrayN1(EV); //status
    plhs[5] = getMexArrayN1(FV); //generation
    plhs[6] = getMexArrayN1(G); //wavelength
    plhs[7] = getMexArrayN1(HV); //parent index
    return;
}