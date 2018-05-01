#include <iostream>
#include <fstream>
#include <cmath>
#include "var_struct.cpp"

using namespace std;

void Roe_flux(double* STL, double* STR, double* n, double* fluxes){

    double FL[5],FR[5],du[5],l[3],gam,gmi,
           rL,uL,vL,wL,unL,qL,pL,rHL,HL,cL,
           rR,uR,vR,wR,unR,qR,pR,rHR,HR,cR,
           di,d1,ui,vi,wi,Hi,af,ucp,c2i,ci,ci1,
           eps,l3,s1,s2,G1,G2,C1,C2;

    gam = 1.4;
    gmi = gam-1.0;

    rL = STL[0]+1e-10;
    
    uL = STL[1]/rL;
    vL = STL[2]/rL;
    wL = STL[3]/rL;

    unL = uL*n[0]+vL*n[1]+wL*n[2];

    qL = sqrt(pow(STL[1],2)+pow(STL[2],2)+pow(STL[3],2))/rL;

    pL = (gam-1.0)*(STL[4]-0.5*rL*qL*qL);

    if ((pL<=0) || (rL<=0)){

        cout<<"Negative density or pressure left cell"<<endl;
    
    }

    rHL = STL[4]+pL;
    HL = rHL/rL;
    cL = sqrt(gam*pL/rL);
    
    FL[0]   = rL*unL;
    FL[1]   = STL[1]*unL+pL*n[0];
    FL[2]   = STL[2]*unL+pL*n[1];
    FL[3]   = STL[3]*unL+pL*n[2];
    FL[4]   = rHL*unL;
    
    ////////////////////////////////////////////////////////////////////////////////////////   
    
    rR = STR[0]+1e-10;

    uR = STR[1]/rR;
    vR = STR[2]/rR;
    wR = STR[3]/rR;

    unR = uR*n[0]+vR*n[1]+wR*n[2];

    qR = sqrt(pow(STR[1],2)+pow(STR[2],2)+pow(STR[3],2))/rR;
    
    pR = (gam-1.0)*(STR[4]-0.5*rR*qR*qR);

    if ((pR<=0) || (rR<=0)){
     
        cout<<"Negative density or pressure right cell"<<endl;
    
    }

    rHR = STR[4]+pR;
    HR = rHR/rR;
    cR = sqrt(gam*pR/rR);
    
    FR[0]   = rR*unR;
    FR[1]   = STR[1]*unR + pR*n[0];
    FR[2]   = STR[2]*unR + pR*n[1];
    FR[3]   = STR[3]*unR + pR*n[2];
    FR[4]   = rHR*unR;
    
    
    ////////////////////////////////////////////////////////////////////////////////////////
    
    for(int j=0; j<5; j++){
        
        du[j]     = STR[j]-STL[j];

    }
    
    di     = sqrt(rR/rL);
    d1     = 1.0/(1.0+di);
    
    ui     = (di*uR+uL)*d1;
    vi     = (di*vR+vL)*d1;
    wi     = (di*wR+wL)*d1;
    Hi     = (di*HR+HL)*d1;
    
    af     = 0.5*(pow(ui,2)+pow(vi,2)+pow(wi,2));
    ucp    = ui*n[0]+vi*n[1]+wi*n[2];
    c2i    = gmi*(Hi-af);
    
    if (c2i<=0){
     
        cout<<"Imaginary Wavespeed"<<endl;
    
    }
    
    ci     = sqrt(c2i);
    ci1    = 1.0/ci;
    
    ////////////////////////////////////////////////////////////////////////////////////////

    l[0]    = ucp+ci;
    l[1]    = ucp-ci;
    l[2]    = ucp;
    
    eps = ci*0.1;
    
    for(int j=0; j<3; j++){
        l[j] = abs(l[j]);
        if (l[j]<eps){
            l[j] = 0.5*(eps + pow(l[j],2)/eps);
        }
    }
    
    l3 = l[2];
    
    s1    = 0.5*(l[0]+l[1]);
    s2    = 0.5*(l[0]-l[1]);
    
    G1    = gmi*(af*du[0]-ui*du[1]-vi*du[2]-wi*du[3]+du[4]);
    G2    = -ucp*du[0]+du[1]*n[0]+du[2]*n[1]+du[3]*n[2];
    
    C1    = G1*(s1-l3)*ci1*ci1+G2*s2*ci1;
    C2    = G1*s2*ci1+G2*(s1-l3);
     
    fluxes[0]    = 0.5*(FL[0]+FR[0])-0.5*(l3*du[0]+C1);
    fluxes[1]    = 0.5*(FL[1]+FR[1])-0.5*(l3*du[1]+C1*ui+C2*n[0]);
    fluxes[2]    = 0.5*(FL[2]+FR[2])-0.5*(l3*du[2]+C1*vi+C2*n[1]);
    fluxes[3]    = 0.5*(FL[3]+FR[3])-0.5*(l3*du[3]+C1*wi+C2*n[2]);
    fluxes[4]    = 0.5*(FL[4]+FR[4])-0.5*(l3*du[4]+C1*Hi+C2*ucp);

}