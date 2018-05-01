#include <iostream>
#include <fstream>
#include "var_struct.cpp"

using namespace std;

void getFlux(Face* face, double* lState, double* rState, BC* bc){

    // Internal face

    if (face->getBCType()==0){

        Roe_flux(lState, rState, face->getNormal(), face->flux);

    }

    // Freestream conditions

    else if (bc[face->getBCType()-1].bctype==1){

        if (face->getLeftCellGID()==-1){

            Roe_flux(bc[face->getBCType()].bcstate, rState, face->getNormal(), face->flux);

        }

        else if (face->getRightCellGID()==-1){

            Roe_flux(lState, bc[face->getBCType()].bcstate, face->getNormal(), face->flux);

        }

    }

    // Inviscid boundary conditions

    else if (bc[face->getBCType()-1].bctype==2){

        double* States;

        if (face->getLeftCellGID()==-1){

            States = rStates;

        }

        else if (face->getRightCellGID()==-1){

            States = lStates;

        }

        double rho  = States[0];
        double rhoU = States[1];
        double rhoV = States[2];
        double rhoW = States[3];
        double rhoE = States[4];
        
        double n1 = face->getNormal()[0];
        double n2 = face->getNormal()[1];
        double n3 = face->getNormal()[2];
        
        double Unm  = (rhoU*n1 + rhoV*n2 + rhoW*n3)/rho;
        double Un   = Unm*n1;
        double Vn   = Unm*n2;
        double Wn   = Unm*n3;
        
        double Ub2  = (rhoU/rho-Un)*(rhoU/rho-Un)+\
                      (rhoV/rho-Vn)*(rhoV/rho-Vn)+\
                      (rhoW/rho-Wn)*(rhoW/rho-Wn);
        
        double pb   = 0.4*(rhoE-0.5*rho*Ub2);

        double flux[5];
        
        flux[0] = 0.0;
        flux[1] = pb*n1;
        flux[2] = pb*n2;
        flux[3] = pb*n3;
        flux[4] = 0.0;

        face->flux = flux;

    }

    // Inflow Conditions

    else if (bc[face->getBCType()-1].bctype==3){

        double* States;

        if (face->getLeftCellGID()==-1){

            States = rStates;

        }

        else if (face->getRightCellGID()==-1){

            States = lStates;

        }

        double rho  = States[0];
        double U    = States[1]/rho;
        double V    = States[2]/rho;
        double W    = States[3]/rho;
        double En   = States[4]/rho;
        double cs   = sqrt(1.4*0.4*(En - 0.5*(U*U+V*V+W*W)));
        
        double n1 = face->getNormal()[0];
        double n2 = face->getNormal()[1];
        double n3 = face->getNormal()[2];
        
        double p0  = bc[face->getBCType()].bcstates[0];
        double ni1 = bc[face->getBCType()].bcstates[1];
        double ni2 = bc[face->getBCType()].bcstates[2];
        double ni3 = bc[face->getBCType()].bcstates[3];
        double T0  = bc[face->getBCType()].bcstates[4];
        
        double Unm = U*n1 + V*n2 + W*n3;
        
        double J  = Unm + (2*cs)/0.4;
        double dn = n1*ni1 + n2*ni2 + n3*ni3;
        
        double aaa = 1.4*287*T0*dn*dn-0.2*J*J;
        double bbb = 4*1.4*287*T0*dn/0.4;
        double ccc = 4*1.4*287*T0/0.16-J*J;
        
        double Mb  = (sqrt(bbb*bbb-4*aaa*ccc)-bbb)/(2*aaa);
        
        double Tb  = T0/(1+0.2*Mb*Mb);
        double pb  = pow(p0*(Tb/T0),1.4/0.4);
        double rb  = pb/(287*Tb);
        double cb  = sqrt(1.4*pb/rb);
        double vb1 = Mb*cb*ni1;
        double vb2 = Mb*cb*ni2;
        double vb3 = Mb*cb*ni3;
        double rEb = pb/0.4 + 0.5*rb*(vb1*vb1+vb2*vb2+vb3*vb3);
        
        double flux[5];
        
        flux[0] = rb*(vb1*n1 + vb2*n2 + vb3*n3);
        flux[1] = rb*vb1*(vb1*n1 + vb2*n2 + vb3*n3)+pb*n1;
        flux[2] = rb*vb2*(vb1*n1 + vb2*n2 + vb3*n3)+pb*n2;
        flux[3] = rb*vb3*(vb1*n1 + vb2*n2 + vb3*n3)+pb*n3;
        flux[4] = (rEb+pb)*(vb1*n1 + vb2*n2 + vb3*n3);

        face->flux = flux;

    }

    // Subsonic Outflow conditions

    else if (bc[face->getBCType()-1].bctype==4){

        double* States;

        if (face->getLeftCellGID()==-1){

            States = rStates;

        }

        else if (face->getRightCellGID()==-1){

            States = lStates;

        }

        double rho  = States[0];
        double U    = States[1]/rho;
        double V    = States[2]/rho;
        double W    = States[3]/rho;
        double En   = States[4]/rho;
        double p    = rho*0.4*(En - 0.5*(U*U+V*V+W*W));
        double cs   = sqrt(1.4*0.4*(En - 0.5*(U*U+V*V+W*W)));
        
        double n1 = face->getNormal()[0];
        double n2 = face->getNormal()[1];
        double n3 = face->getNormal()[2];

        double pb  = bc[face->getBCType()].bcstates[0];
        double rb  = pow(pb*pow(rho,1.4)/p,1/1.4);
        double cb  = sqrt(1.4*pb/rb);

        double Unm = U*n1 + V*n2 + W*n3;
        double J  = Unm + (2*cs)/0.4;

        double Unbm = J - 2*cb/0.4;

        double vb1 = U-Unm*n1+Unbm*n1;
        double vb2 = V-Unm*n2+Unbm*n2;
        double vb3 = W-Unm*n3+Unbm*n3;

        double rEb = pb/0.4 + 0.5*rb*(vb1*vb1+vb2*vb2+vb3*vb3);

        double flux[5];

        flux[0] = rb*(vb1*n1 + vb2*n2 + vb3*n3);
        flux[1] = rb*vb1*(vb1*n1 + vb2*n2 + vb3*n3)+pb*n1;
        flux[2] = rb*vb2*(vb1*n1 + vb2*n2 + vb3*n3)+pb*n2;
        flux[3] = rb*vb3*(vb1*n1 + vb2*n2 + vb3*n3)+pb*n3;
        flux[4] = (rEb+pb)*(vb1*n1 + vb2*n2 + vb3*n3);

        face->flux = flux;

    }

}