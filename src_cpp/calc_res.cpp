#include <iostream>
#include <cmath>
#include "flux.cpp"

using namespace std;



void bc_flux(int nblock, block* mesh){

    for(int block_id=0; block_id<nblock; block_id++){
    
        for(int bc=0; bc<6; bc++){
        
            // Periodic boundary condition

            if (mesh[block_id].bc_type[bc]<=0){
            

                if (bc==0){
                
                    for(int j=0; j<mesh[block_id].ny-1; j++){
                    
                        for(int k=0; k<mesh[block_id].nz-1; k++){
                        
                            Roe_flux(mesh[ -mesh[block_id].bc_type[bc] ].states[ mesh[ -mesh[block_id].bc_type[bc] ].nx-1 ][j][k],\
                                     mesh[ block_id ].states[0][j][k],\
                                     mesh[ block_id ].xn[0][j][k],\
                                     mesh[ block_id ].xfluxes[0][j][k]);
                        
                        }
                    
                    }
                
                }
            
                if (bc==1){
                
                    for(int j=0; j<mesh[block_id].ny-1; j++){
                    
                        for(int k=0; k<mesh[block_id].nz-1; k++){
                        
                            Roe_flux(mesh[ block_id ].states[ mesh[block_id].nx-1 ][j][k],\
                                     mesh[ -mesh[block_id].bc_type[bc] ].states[0][j][k],\
                                     mesh[ block_id ].xn[ mesh[block_id].nx ][j][k],\
                                     mesh[ block_id ].xfluxes[ mesh[block_id].nx ][j][k]);
                        
                        }
                    
                    }
                
                }
            
                if (bc==2){
                
                    for(int i=0; i<mesh[block_id].nx-1; i++){
                    
                        for(int k=0; k<mesh[block_id].nz-1; k++){
                        
                            Roe_flux(mesh[ -mesh[block_id].bc_type[bc] ].states[i][ mesh[ -mesh[block_id].bc_type[bc] ].ny-1 ][k],\
                                     mesh[ block_id ].states[i][0][k],\
                                     mesh[ block_id ].yn[i][0][k],\
                                     mesh[ block_id ].yfluxes[i][0][k]);
                        
                        }
                    
                    }
                
                }
            
                if (bc==3){
                
                    for(int i=0; i<mesh[block_id].nx-1; i++){
                    
                        for(int k=0; k<mesh[block_id].nz-1; k++){
                        
                            Roe_flux(mesh[ block_id ].states[i][ mesh[block_id].ny-1 ][k],\
                                     mesh[ -mesh[block_id].bc_type[bc] ].states[i][0][k],\
                                     mesh[ block_id ].yn[i][ mesh[block_id].ny ][k],\
                                     mesh[ block_id ].yfluxes[i][ mesh[block_id].ny ][k]);
                        
                        }
                    
                    }
                
                }
            
                if (bc==4){
                
                    for(int i=0; i<mesh[block_id].nx-1; i++){
                    
                        for(int j=0; j<mesh[block_id].ny-1; j++){
                        
                            Roe_flux(mesh[ -mesh[block_id].bc_type[bc] ].states[i][j][ mesh[ -mesh[block_id].bc_type[bc] ].nz-1 ],\
                                     mesh[ block_id ].states[i][j][0],\
                                     mesh[ block_id ].zn[i][j][0],\
                                     mesh[ block_id ].zfluxes[i][j][0]);
                        
                        }
                    
                    }
                
                }
            
                if (bc==5){
                
                    for(int i=0; i<mesh[block_id].nx-1; i++){
                    
                        for(int j=0; j<mesh[block_id].ny-1; j++){
                        
                            Roe_flux(mesh[ block_id ].states[i][j][ mesh[block_id].nz-1 ],\
                                     mesh[ -mesh[block_id].bc_type[bc] ].states[i][j][0],\
                                     mesh[ block_id ].zn[i][j][ mesh[block_id].nz ],\
                                     mesh[ block_id ].zfluxes[i][j][ mesh[block_id].nz ]);
                        
                        }
                    
                    }
                
                }

            
            }

            // Inviscid wall

            if(mesh[block_id].bc_type[bc]==1){

                if (bc==0){
                
                    for(int j=0; j<mesh[block_id].ny-1; j++){
                    
                        for(int k=0; k<mesh[block_id].nz-1; k++){
                        
                            int c = 0;
                            int f = 0;

                            double rho  = mesh[block_id].states[c][j][k][0];
                            double rhoU = mesh[block_id].states[c][j][k][1];
                            double rhoV = mesh[block_id].states[c][j][k][2];
                            double rhoW = mesh[block_id].states[c][j][k][3];
                            double rhoE = mesh[block_id].states[c][j][k][4];

                            double n1   = mesh[block_id].xn[f][j][k][0];
                            double n2   = mesh[block_id].xn[f][j][k][1];
                            double n3   = mesh[block_id].xn[f][j][k][2];

                            double Unm  = (rhoU*n1 + rhoV*n2 + rhoW*n3)/rho;
                            double Un   = Unm*n1;
                            double Vn   = Unm*n2;
                            double Wn   = Unm*n3;
                            double Ub2  = (rhoU/rho-Un)*(rhoU/rho-Un)+\
                                          (rhoV/rho-Vn)*(rhoV/rho-Vn)+\
                                          (rhoW/rho-Wn)*(rhoW/rho-Wn);
                            
                            double pb   = 0.4*(rhoE-0.5*rho*Ub2);

                            mesh[block_id].xfluxes[f][j][k][0] = 0.0;
                            mesh[block_id].xfluxes[f][j][k][1] = pb*n1;
                            mesh[block_id].xfluxes[f][j][k][2] = pb*n2;
                            mesh[block_id].xfluxes[f][j][k][3] = pb*n3;
                            mesh[block_id].xfluxes[f][j][k][4] = 0.0;

                        }
                    
                    }
                
                }
            
                if (bc==1){
                
                    for(int j=0; j<mesh[block_id].ny-1; j++){
                    
                        for(int k=0; k<mesh[block_id].nz-1; k++){
                        
                            int c = mesh[block_id].nx-1;
                            int f = c+1;

                            double rho  = mesh[block_id].states[c][j][k][0];
                            double rhoU = mesh[block_id].states[c][j][k][1];
                            double rhoV = mesh[block_id].states[c][j][k][2];
                            double rhoW = mesh[block_id].states[c][j][k][3];
                            double rhoE = mesh[block_id].states[c][j][k][4];

                            double n1   = mesh[block_id].xn[f][j][k][0];
                            double n2   = mesh[block_id].xn[f][j][k][1];
                            double n3   = mesh[block_id].xn[f][j][k][2];

                            double Unm  = (rhoU*n1 + rhoV*n2 + rhoW*n3)/rho;
                            double Un   = Unm*n1;
                            double Vn   = Unm*n2;
                            double Wn   = Unm*n3;
                            double Ub2  = (rhoU/rho-Un)*(rhoU/rho-Un)+\
                                          (rhoV/rho-Vn)*(rhoV/rho-Vn)+\
                                          (rhoW/rho-Wn)*(rhoW/rho-Wn);
                            
                            double pb   = 0.4*(rhoE-0.5*rho*Ub2);

                            mesh[block_id].xfluxes[f][j][k][0] = 0.0;
                            mesh[block_id].xfluxes[f][j][k][1] = pb*n1;
                            mesh[block_id].xfluxes[f][j][k][2] = pb*n2;
                            mesh[block_id].xfluxes[f][j][k][3] = pb*n3;
                            mesh[block_id].xfluxes[f][j][k][4] = 0.0;
                        
                        }
                    
                    }
                
                }
            
                if (bc==2){
                
                    for(int i=0; i<mesh[block_id].nx-1; i++){
                    
                        for(int k=0; k<mesh[block_id].nz-1; k++){
                        
                            int c = 0;
                            int f = 0;

                            double rho  = mesh[block_id].states[i][c][k][0];
                            double rhoU = mesh[block_id].states[i][c][k][1];
                            double rhoV = mesh[block_id].states[i][c][k][2];
                            double rhoW = mesh[block_id].states[i][c][k][3];
                            double rhoE = mesh[block_id].states[i][c][k][4];

                            double n1   = mesh[block_id].yn[i][f][k][0];
                            double n2   = mesh[block_id].yn[i][f][k][1];
                            double n3   = mesh[block_id].yn[i][f][k][2];

                            double Unm  = (rhoU*n1 + rhoV*n2 + rhoW*n3)/rho;
                            double Un   = Unm*n1;
                            double Vn   = Unm*n2;
                            double Wn   = Unm*n3;
                            double Ub2  = (rhoU/rho-Un)*(rhoU/rho-Un)+\
                                          (rhoV/rho-Vn)*(rhoV/rho-Vn)+\
                                          (rhoW/rho-Wn)*(rhoW/rho-Wn);
                            
                            double pb   = 0.4*(rhoE-0.5*rho*Ub2);

                            mesh[block_id].yfluxes[i][f][k][0] = 0.0;
                            mesh[block_id].yfluxes[i][f][k][1] = pb*n1;
                            mesh[block_id].yfluxes[i][f][k][2] = pb*n2;
                            mesh[block_id].yfluxes[i][f][k][3] = pb*n3;
                            mesh[block_id].yfluxes[i][f][k][4] = 0.0;
                        
                        }
                    
                    }
                
                }
            
                if (bc==3){
                
                    for(int i=0; i<mesh[block_id].nx-1; i++){
                    
                        for(int k=0; k<mesh[block_id].nz-1; k++){
                        
                            int c = mesh[block_id].ny-1;
                            int f = c+1;

                            double rho  = mesh[block_id].states[i][c][k][0];
                            double rhoU = mesh[block_id].states[i][c][k][1];
                            double rhoV = mesh[block_id].states[i][c][k][2];
                            double rhoW = mesh[block_id].states[i][c][k][3];
                            double rhoE = mesh[block_id].states[i][c][k][4];

                            double n1   = mesh[block_id].yn[i][f][k][0];
                            double n2   = mesh[block_id].yn[i][f][k][1];
                            double n3   = mesh[block_id].yn[i][f][k][2];

                            double Unm  = (rhoU*n1 + rhoV*n2 + rhoW*n3)/rho;
                            double Un   = Unm*n1;
                            double Vn   = Unm*n2;
                            double Wn   = Unm*n3;
                            double Ub2  = (rhoU/rho-Un)*(rhoU/rho-Un)+\
                                          (rhoV/rho-Vn)*(rhoV/rho-Vn)+\
                                          (rhoW/rho-Wn)*(rhoW/rho-Wn);
                            
                            double pb   = 0.4*(rhoE-0.5*rho*Ub2);

                            mesh[block_id].yfluxes[i][f][k][0] = 0.0;
                            mesh[block_id].yfluxes[i][f][k][1] = pb*n1;
                            mesh[block_id].yfluxes[i][f][k][2] = pb*n2;
                            mesh[block_id].yfluxes[i][f][k][3] = pb*n3;
                            mesh[block_id].yfluxes[i][f][k][4] = 0.0;
                        
                        }
                    
                    }
                
                }
            
                if (bc==4){
                
                    for(int i=0; i<mesh[block_id].nx-1; i++){
                    
                        for(int j=0; j<mesh[block_id].ny-1; j++){
                        
                            int c = 0;
                            int f = 0;

                            double rho  = mesh[block_id].states[i][j][c][0];
                            double rhoU = mesh[block_id].states[i][j][c][1];
                            double rhoV = mesh[block_id].states[i][j][c][2];
                            double rhoW = mesh[block_id].states[i][j][c][3];
                            double rhoE = mesh[block_id].states[i][j][c][4];

                            double n1   = mesh[block_id].zn[i][j][f][0];
                            double n2   = mesh[block_id].zn[i][j][f][1];
                            double n3   = mesh[block_id].zn[i][j][f][2];

                            double Unm  = (rhoU*n1 + rhoV*n2 + rhoW*n3)/rho;
                            double Un   = Unm*n1;
                            double Vn   = Unm*n2;
                            double Wn   = Unm*n3;
                            double Ub2  = (rhoU/rho-Un)*(rhoU/rho-Un)+\
                                          (rhoV/rho-Vn)*(rhoV/rho-Vn)+\
                                          (rhoW/rho-Wn)*(rhoW/rho-Wn);
                            
                            double pb   = 0.4*(rhoE-0.5*rho*Ub2);

                            mesh[block_id].zfluxes[i][j][f][0] = 0.0;
                            mesh[block_id].zfluxes[i][j][f][1] = pb*n1;
                            mesh[block_id].zfluxes[i][j][f][2] = pb*n2;
                            mesh[block_id].zfluxes[i][j][f][3] = pb*n3;
                            mesh[block_id].zfluxes[i][j][f][4] = 0.0;
                        
                        }
                    
                    }
                
                }
            
                if (bc==5){
                
                    for(int i=0; i<mesh[block_id].nx-1; i++){
                    
                        for(int j=0; j<mesh[block_id].ny-1; j++){
                        
                            int c = mesh[block_id].nz-1;
                            int f = c+1;

                            double rho  = mesh[block_id].states[i][j][c][0];
                            double rhoU = mesh[block_id].states[i][j][c][1];
                            double rhoV = mesh[block_id].states[i][j][c][2];
                            double rhoW = mesh[block_id].states[i][j][c][3];
                            double rhoE = mesh[block_id].states[i][j][c][4];

                            double n1   = mesh[block_id].zn[i][j][f][0];
                            double n2   = mesh[block_id].zn[i][j][f][1];
                            double n3   = mesh[block_id].zn[i][j][f][2];

                            double Unm  = (rhoU*n1 + rhoV*n2 + rhoW*n3)/rho;
                            double Un   = Unm*n1;
                            double Vn   = Unm*n2;
                            double Wn   = Unm*n3;
                            double Ub2  = (rhoU/rho-Un)*(rhoU/rho-Un)+\
                                          (rhoV/rho-Vn)*(rhoV/rho-Vn)+\
                                          (rhoW/rho-Wn)*(rhoW/rho-Wn);
                            
                            double pb   = 0.4*(rhoE-0.5*rho*Ub2);

                            mesh[block_id].zfluxes[i][j][f][0] = 0.0;
                            mesh[block_id].zfluxes[i][j][f][1] = pb*n1;
                            mesh[block_id].zfluxes[i][j][f][2] = pb*n2;
                            mesh[block_id].zfluxes[i][j][f][3] = pb*n3;
                            mesh[block_id].zfluxes[i][j][f][4] = 0.0;
                        
                        }
                    
                    }
                
                }

            }

            // Far-field condition

            if(mesh[block_id].bc_type[bc]==2){

                if (bc==0){
                
                    for(int j=0; j<mesh[block_id].ny-1; j++){
                    
                        for(int k=0; k<mesh[block_id].nz-1; k++){
                        
                            int c = 0;
                            int f = 0;

                            Roe_flux(mesh[block_id].bc_states[bc],\
                                     mesh[block_id].states[c][j][k],\
                                     mesh[block_id].xn[f][j][k],\
                                     mesh[block_id].xfluxes[f][j][k]);            
                            
                        }
                    
                    }
                
                }
            
                if (bc==1){
                
                    for(int j=0; j<mesh[block_id].ny-1; j++){
                    
                        for(int k=0; k<mesh[block_id].nz-1; k++){
                        
                            int c = mesh[block_id].nx-1;
                            int f = c+1;

                            Roe_flux(mesh[block_id].states[c][j][k],\
                                     mesh[block_id].bc_states[bc],\
                                     mesh[block_id].xn[f][j][k],\
                                     mesh[block_id].xfluxes[f][j][k]);

                        }
                    
                    }
                
                }
            
                if (bc==2){
                
                    for(int i=0; i<mesh[block_id].nx-1; i++){
                    
                        for(int k=0; k<mesh[block_id].nz-1; k++){
                        
                            int c = 0;
                            int f = 0;

                            Roe_flux(mesh[block_id].bc_states[bc],\
                                     mesh[block_id].states[i][c][k],\
                                     mesh[block_id].yn[i][f][k],\
                                     mesh[block_id].yfluxes[i][f][k]);
                        
                        }
                    
                    }
                
                }
            
                if (bc==3){
                
                    for(int i=0; i<mesh[block_id].nx-1; i++){
                    
                        for(int k=0; k<mesh[block_id].nz-1; k++){
                        
                            int c = mesh[block_id].ny-1;
                            int f = c+1;

                            Roe_flux(mesh[block_id].states[i][c][k],\
                                     mesh[block_id].bc_states[bc],\
                                     mesh[block_id].yn[i][f][k],\
                                     mesh[block_id].yfluxes[i][f][k]);
                        
                        }
                    
                    }
                
                }
            
                if (bc==4){
                
                    for(int i=0; i<mesh[block_id].nx-1; i++){
                    
                        for(int j=0; j<mesh[block_id].ny-1; j++){
                        
                            int c = 0;
                            int f = 0;

                            Roe_flux(mesh[block_id].bc_states[bc],\
                                     mesh[block_id].states[i][j][c],\
                                     mesh[block_id].zn[i][j][f],\
                                     mesh[block_id].zfluxes[i][j][f]);
                        
                        }
                    
                    }
                
                }
            
                if (bc==5){
                
                    for(int i=0; i<mesh[block_id].nx-1; i++){
                    
                        for(int j=0; j<mesh[block_id].ny-1; j++){
                        
                            int c = mesh[block_id].nz-1;
                            int f = c+1;

                            Roe_flux(mesh[block_id].states[i][j][c],\
                                     mesh[block_id].bc_states[bc],\
                                     mesh[block_id].zn[i][j][f],\
                                     mesh[block_id].zfluxes[i][j][f]);
                        
                        }
                    
                    }
                
                }

            }

            // Inflow condition

            if(mesh[block_id].bc_type[bc]==3){

                if (bc==0){
                
                    for(int j=0; j<mesh[block_id].ny-1; j++){
                    
                        for(int k=0; k<mesh[block_id].nz-1; k++){
                        
                            int c = 0;
                            int f = 0;

                            double rho = mesh[block_id].states[c][j][k][0];
                            double U = mesh[block_id].states[c][j][k][1]/rho;
                            double V = mesh[block_id].states[c][j][k][2]/rho;
                            double W = mesh[block_id].states[c][j][k][3]/rho;
                            double En = mesh[block_id].states[c][j][k][4]/rho;
                            double cs = sqrt(1.4*0.4*(En - 0.5*(U*U+V*V+W*W)));

                            double n1 = mesh[block_id].xn[f][j][k][0];
                            double n2 = mesh[block_id].xn[f][j][k][1];
                            double n3 = mesh[block_id].xn[f][j][k][2];

                            double p0  = mesh[block_id].bc_states[bc][0];
                            double ni1 = mesh[block_id].bc_states[bc][1];
                            double ni2 = mesh[block_id].bc_states[bc][2];
                            double ni3 = mesh[block_id].bc_states[bc][3];
                            double T0  = mesh[block_id].bc_states[bc][4];

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

                            mesh[block_id].xfluxes[f][j][k][0] = flux[0];
                            mesh[block_id].xfluxes[f][j][k][1] = flux[1];
                            mesh[block_id].xfluxes[f][j][k][2] = flux[2];
                            mesh[block_id].xfluxes[f][j][k][3] = flux[3];
                            mesh[block_id].xfluxes[f][j][k][4] = flux[4];

                        }
                    
                    }
                
                }
            
                if (bc==1){
                
                    for(int j=0; j<mesh[block_id].ny-1; j++){
                    
                        for(int k=0; k<mesh[block_id].nz-1; k++){
                        
                            int c = mesh[block_id].nx-1;
                            int f = c+1;

                            double rho = mesh[block_id].states[c][j][k][0];
                            double U = mesh[block_id].states[c][j][k][1]/rho;
                            double V = mesh[block_id].states[c][j][k][2]/rho;
                            double W = mesh[block_id].states[c][j][k][3]/rho;
                            double En = mesh[block_id].states[c][j][k][4]/rho;
                            double cs = sqrt(1.4*0.4*(En - 0.5*(U*U+V*V+W*W)));

                            double n1 = mesh[block_id].xn[f][j][k][0];
                            double n2 = mesh[block_id].xn[f][j][k][1];
                            double n3 = mesh[block_id].xn[f][j][k][2];

                            double p0  = mesh[block_id].bc_states[bc][0];
                            double ni1 = mesh[block_id].bc_states[bc][1];
                            double ni2 = mesh[block_id].bc_states[bc][2];
                            double ni3 = mesh[block_id].bc_states[bc][3];
                            double T0  = mesh[block_id].bc_states[bc][4];

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

                            mesh[block_id].xfluxes[f][j][k][0] = flux[0];
                            mesh[block_id].xfluxes[f][j][k][1] = flux[1];
                            mesh[block_id].xfluxes[f][j][k][2] = flux[2];
                            mesh[block_id].xfluxes[f][j][k][3] = flux[3];
                            mesh[block_id].xfluxes[f][j][k][4] = flux[4];

                        }
                    
                    }
                
                }
            
                if (bc==2){
                
                    for(int i=0; i<mesh[block_id].nx-1; i++){
                    
                        for(int k=0; k<mesh[block_id].nz-1; k++){
                        
                            int c = 0;
                            int f = 0;

                            double rho = mesh[block_id].states[i][c][k][0];
                            double U   = mesh[block_id].states[i][c][k][1]/rho;
                            double V   = mesh[block_id].states[i][c][k][2]/rho;
                            double W   = mesh[block_id].states[i][c][k][3]/rho;
                            double En  = mesh[block_id].states[i][c][k][4]/rho;
                            double cs  = sqrt(1.4*0.4*(En - 0.5*(U*U+V*V+W*W)));

                            double n1 = mesh[block_id].yn[i][f][k][0];
                            double n2 = mesh[block_id].yn[i][f][k][1];
                            double n3 = mesh[block_id].yn[i][f][k][2];

                            double p0  = mesh[block_id].bc_states[bc][0];
                            double ni1 = mesh[block_id].bc_states[bc][1];
                            double ni2 = mesh[block_id].bc_states[bc][2];
                            double ni3 = mesh[block_id].bc_states[bc][3];
                            double T0  = mesh[block_id].bc_states[bc][4];

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

                            mesh[block_id].yfluxes[i][f][k][0] = flux[0];
                            mesh[block_id].yfluxes[i][f][k][1] = flux[1];
                            mesh[block_id].yfluxes[i][f][k][2] = flux[2];
                            mesh[block_id].yfluxes[i][f][k][3] = flux[3];
                            mesh[block_id].yfluxes[i][f][k][4] = flux[4];

                        }
                    
                    }
                
                }
            
                if (bc==3){
                
                    for(int i=0; i<mesh[block_id].nx-1; i++){
                    
                        for(int k=0; k<mesh[block_id].nz-1; k++){
                        
                            int c = mesh[block_id].ny-1;
                            int f = c+1;

                            double rho = mesh[block_id].states[i][c][k][0];
                            double U   = mesh[block_id].states[i][c][k][1]/rho;
                            double V   = mesh[block_id].states[i][c][k][2]/rho;
                            double W   = mesh[block_id].states[i][c][k][3]/rho;
                            double En  = mesh[block_id].states[i][c][k][4]/rho;
                            double cs  = sqrt(1.4*0.4*(En - 0.5*(U*U+V*V+W*W)));

                            double n1 = mesh[block_id].yn[i][f][k][0];
                            double n2 = mesh[block_id].yn[i][f][k][1];
                            double n3 = mesh[block_id].yn[i][f][k][2];

                            double p0  = mesh[block_id].bc_states[bc][0];
                            double ni1 = mesh[block_id].bc_states[bc][1];
                            double ni2 = mesh[block_id].bc_states[bc][2];
                            double ni3 = mesh[block_id].bc_states[bc][3];
                            double T0  = mesh[block_id].bc_states[bc][4];

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

                            mesh[block_id].yfluxes[i][f][k][0] = flux[0];
                            mesh[block_id].yfluxes[i][f][k][1] = flux[1];
                            mesh[block_id].yfluxes[i][f][k][2] = flux[2];
                            mesh[block_id].yfluxes[i][f][k][3] = flux[3];
                            mesh[block_id].yfluxes[i][f][k][4] = flux[4];

                        }
                    
                    }
                
                }
            
                if (bc==4){
                
                    for(int i=0; i<mesh[block_id].nx-1; i++){
                    
                        for(int j=0; j<mesh[block_id].ny-1; j++){
                        
                            int c = 0;
                            int f = 0;

                            double rho = mesh[block_id].states[i][j][c][0];
                            double U   = mesh[block_id].states[i][j][c][1]/rho;
                            double V   = mesh[block_id].states[i][j][c][2]/rho;
                            double W   = mesh[block_id].states[i][j][c][3]/rho;
                            double En  = mesh[block_id].states[i][j][c][4]/rho;
                            double cs  = sqrt(1.4*0.4*(En - 0.5*(U*U+V*V+W*W)));

                            double n1 = mesh[block_id].zn[i][j][f][0];
                            double n2 = mesh[block_id].zn[i][j][f][1];
                            double n3 = mesh[block_id].zn[i][j][f][2];

                            double p0  = mesh[block_id].bc_states[bc][0];
                            double ni1 = mesh[block_id].bc_states[bc][1];
                            double ni2 = mesh[block_id].bc_states[bc][2];
                            double ni3 = mesh[block_id].bc_states[bc][3];
                            double T0  = mesh[block_id].bc_states[bc][4];

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

                            mesh[block_id].zfluxes[i][j][f][0] = flux[0];
                            mesh[block_id].zfluxes[i][j][f][1] = flux[1];
                            mesh[block_id].zfluxes[i][j][f][2] = flux[2];
                            mesh[block_id].zfluxes[i][j][f][3] = flux[3];
                            mesh[block_id].zfluxes[i][j][f][4] = flux[4];

                        }
                    
                    }
                
                }
            
                if (bc==5){
                
                    for(int i=0; i<mesh[block_id].nx-1; i++){
                    
                        for(int j=0; j<mesh[block_id].ny-1; j++){
                        
                            int c = mesh[block_id].nz-1;
                            int f = c+1;

                            double rho = mesh[block_id].states[i][j][c][0];
                            double U   = mesh[block_id].states[i][j][c][1]/rho;
                            double V   = mesh[block_id].states[i][j][c][2]/rho;
                            double W   = mesh[block_id].states[i][j][c][3]/rho;
                            double En  = mesh[block_id].states[i][j][c][4]/rho;
                            double cs  = sqrt(1.4*0.4*(En - 0.5*(U*U+V*V+W*W)));

                            double n1 = mesh[block_id].zn[i][j][f][0];
                            double n2 = mesh[block_id].zn[i][j][f][1];
                            double n3 = mesh[block_id].zn[i][j][f][2];

                            double p0  = mesh[block_id].bc_states[bc][0];
                            double ni1 = mesh[block_id].bc_states[bc][1];
                            double ni2 = mesh[block_id].bc_states[bc][2];
                            double ni3 = mesh[block_id].bc_states[bc][3];
                            double T0  = mesh[block_id].bc_states[bc][4];

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

                            mesh[block_id].zfluxes[i][j][f][0] = flux[0];
                            mesh[block_id].zfluxes[i][j][f][1] = flux[1];
                            mesh[block_id].zfluxes[i][j][f][2] = flux[2];
                            mesh[block_id].zfluxes[i][j][f][3] = flux[3];
                            mesh[block_id].zfluxes[i][j][f][4] = flux[4];

                        }
                    
                    }
                
                }

            }

            // Subsonic outflow

            if(mesh[block_id].bc_type[bc]==4){

                if (bc==0){
                
                    for(int j=0; j<mesh[block_id].ny-1; j++){
                    
                        for(int k=0; k<mesh[block_id].nz-1; k++){
                        
                            int c = 0;
                            int f = 0;

                            double rho = mesh[block_id].states[c][j][k][0];
                            double U   = mesh[block_id].states[c][j][k][1]/rho;
                            double V   = mesh[block_id].states[c][j][k][2]/rho;
                            double W   = mesh[block_id].states[c][j][k][3]/rho;
                            double En  = mesh[block_id].states[c][j][k][4]/rho;
                            double p  = rho*0.4*(En - 0.5*(U*U+V*V+W*W));
                            double cs = sqrt(1.4*p/rho);

                            double n1 = mesh[block_id].xn[f][j][k][0];
                            double n2 = mesh[block_id].xn[f][j][k][1];
                            double n3 = mesh[block_id].xn[f][j][k][2];

                            double pb  = mesh[block_id].bc_states[bc][0];
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

                            mesh[block_id].xfluxes[f][j][k][0] = flux[0];
                            mesh[block_id].xfluxes[f][j][k][1] = flux[1];
                            mesh[block_id].xfluxes[f][j][k][2] = flux[2];
                            mesh[block_id].xfluxes[f][j][k][3] = flux[3];
                            mesh[block_id].xfluxes[f][j][k][4] = flux[4];

                        }
                    
                    }
                
                }
            
                if (bc==1){
                
                    for(int j=0; j<mesh[block_id].ny-1; j++){
                    
                        for(int k=0; k<mesh[block_id].nz-1; k++){
                        
                            int c = mesh[block_id].nx-1;
                            int f = c+1;

                            double rho = mesh[block_id].states[c][j][k][0];
                            double U   = mesh[block_id].states[c][j][k][1]/rho;
                            double V   = mesh[block_id].states[c][j][k][2]/rho;
                            double W   = mesh[block_id].states[c][j][k][3]/rho;
                            double En  = mesh[block_id].states[c][j][k][4]/rho;
                            double p  = rho*0.4*(En - 0.5*(U*U+V*V+W*W));
                            double cs = sqrt(1.4*p/rho);

                            double n1 = mesh[block_id].xn[f][j][k][0];
                            double n2 = mesh[block_id].xn[f][j][k][1];
                            double n3 = mesh[block_id].xn[f][j][k][2];

                            double pb  = mesh[block_id].bc_states[bc][0];
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

                            mesh[block_id].xfluxes[f][j][k][0] = flux[0];
                            mesh[block_id].xfluxes[f][j][k][1] = flux[1];
                            mesh[block_id].xfluxes[f][j][k][2] = flux[2];
                            mesh[block_id].xfluxes[f][j][k][3] = flux[3];
                            mesh[block_id].xfluxes[f][j][k][4] = flux[4];

                        }
                    
                    }
                
                }
            
                if (bc==2){
                
                    for(int i=0; i<mesh[block_id].nx-1; i++){
                    
                        for(int k=0; k<mesh[block_id].nz-1; k++){
                        
                            int c = 0;
                            int f = 0;

                            double rho = mesh[block_id].states[i][c][k][0];
                            double U   = mesh[block_id].states[i][c][k][1]/rho;
                            double V   = mesh[block_id].states[i][c][k][2]/rho;
                            double W   = mesh[block_id].states[i][c][k][3]/rho;
                            double En  = mesh[block_id].states[i][c][k][4]/rho;
                            double p  = rho*0.4*(En - 0.5*(U*U+V*V+W*W));
                            double cs = sqrt(1.4*p/rho);

                            double n1 = mesh[block_id].yn[i][f][k][0];
                            double n2 = mesh[block_id].yn[i][f][k][1];
                            double n3 = mesh[block_id].yn[i][f][k][2];

                            double pb  = mesh[block_id].bc_states[bc][0];
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

                            mesh[block_id].yfluxes[i][f][k][0] = flux[0];
                            mesh[block_id].yfluxes[i][f][k][1] = flux[1];
                            mesh[block_id].yfluxes[i][f][k][2] = flux[2];
                            mesh[block_id].yfluxes[i][f][k][3] = flux[3];
                            mesh[block_id].yfluxes[i][f][k][4] = flux[4];

                        }
                    
                    }
                
                }
            
                if (bc==3){
                
                    for(int i=0; i<mesh[block_id].nx-1; i++){
                    
                        for(int k=0; k<mesh[block_id].nz-1; k++){
                        
                            int c = mesh[block_id].ny-1;
                            int f = c+1;

                            double rho = mesh[block_id].states[i][c][k][0];
                            double U   = mesh[block_id].states[i][c][k][1]/rho;
                            double V   = mesh[block_id].states[i][c][k][2]/rho;
                            double W   = mesh[block_id].states[i][c][k][3]/rho;
                            double En  = mesh[block_id].states[i][c][k][4]/rho;
                            double p  = rho*0.4*(En - 0.5*(U*U+V*V+W*W));
                            double cs = sqrt(1.4*p/rho);

                            double n1 = mesh[block_id].yn[i][f][k][0];
                            double n2 = mesh[block_id].yn[i][f][k][1];
                            double n3 = mesh[block_id].yn[i][f][k][2];

                            double pb  = mesh[block_id].bc_states[bc][0];
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

                            mesh[block_id].yfluxes[i][f][k][0] = flux[0];
                            mesh[block_id].yfluxes[i][f][k][1] = flux[1];
                            mesh[block_id].yfluxes[i][f][k][2] = flux[2];
                            mesh[block_id].yfluxes[i][f][k][3] = flux[3];
                            mesh[block_id].yfluxes[i][f][k][4] = flux[4];

                        }
                    
                    }
                
                }
            
                if (bc==4){
                
                    for(int i=0; i<mesh[block_id].nx-1; i++){
                    
                        for(int j=0; j<mesh[block_id].ny-1; j++){
                        
                            int c = 0;
                            int f = 0;

                            double rho = mesh[block_id].states[i][j][c][0];
                            double U   = mesh[block_id].states[i][j][c][1]/rho;
                            double V   = mesh[block_id].states[i][j][c][2]/rho;
                            double W   = mesh[block_id].states[i][j][c][3]/rho;
                            double En  = mesh[block_id].states[i][j][c][4]/rho;
                            double p  = rho*0.4*(En - 0.5*(U*U+V*V+W*W));
                            double cs = sqrt(1.4*p/rho);

                            double n1 = mesh[block_id].zn[i][j][f][0];
                            double n2 = mesh[block_id].zn[i][j][f][1];
                            double n3 = mesh[block_id].zn[i][j][f][2];

                            double pb  = mesh[block_id].bc_states[bc][0];
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

                            mesh[block_id].zfluxes[i][j][f][0] = flux[0];
                            mesh[block_id].zfluxes[i][j][f][1] = flux[1];
                            mesh[block_id].zfluxes[i][j][f][2] = flux[2];
                            mesh[block_id].zfluxes[i][j][f][3] = flux[3];
                            mesh[block_id].zfluxes[i][j][f][4] = flux[4];

                        }
                    
                    }
                
                }
            
                if (bc==5){
                
                    for(int i=0; i<mesh[block_id].nx-1; i++){
                    
                        for(int j=0; j<mesh[block_id].ny-1; j++){
                        
                            int c = mesh[block_id].nz-1;
                            int f = c+1;

                            double rho = mesh[block_id].states[i][j][c][0];
                            double U   = mesh[block_id].states[i][j][c][1]/rho;
                            double V   = mesh[block_id].states[i][j][c][2]/rho;
                            double W   = mesh[block_id].states[i][j][c][3]/rho;
                            double En  = mesh[block_id].states[i][j][c][4]/rho;
                            double p  = rho*0.4*(En - 0.5*(U*U+V*V+W*W));
                            double cs = sqrt(1.4*p/rho);

                            double n1 = mesh[block_id].zn[i][j][f][0];
                            double n2 = mesh[block_id].zn[i][j][f][1];
                            double n3 = mesh[block_id].zn[i][j][f][2];

                            double pb  = mesh[block_id].bc_states[bc][0];
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

                            mesh[block_id].zfluxes[i][j][f][0] = flux[0];
                            mesh[block_id].zfluxes[i][j][f][1] = flux[1];
                            mesh[block_id].zfluxes[i][j][f][2] = flux[2];
                            mesh[block_id].zfluxes[i][j][f][3] = flux[3];
                            mesh[block_id].zfluxes[i][j][f][4] = flux[4];

                        }
                    
                    }
                
                }

            }

            // Supersonic outflow

            if(mesh[block_id].bc_type[bc]==5){

                if (bc==0){
                
                    for(int j=0; j<mesh[block_id].ny-1; j++){
                    
                        for(int k=0; k<mesh[block_id].nz-1; k++){
                        
                            int c = 0;
                            int f = 0;

                            double rb  = mesh[block_id].states[c][j][k][0];
                            double vb1 = mesh[block_id].states[c][j][k][1]/rb;
                            double vb2 = mesh[block_id].states[c][j][k][2]/rb;
                            double vb3 = mesh[block_id].states[c][j][k][3]/rb;
                            double rEb = mesh[block_id].states[c][j][k][4];
                            double pb  = 0.4*(rEb - 0.5*rb*(vb1*vb1+vb2*vb2+vb3*vb3));
                            
                            double n1 = mesh[block_id].xn[f][j][k][0];
                            double n2 = mesh[block_id].xn[f][j][k][1];
                            double n3 = mesh[block_id].xn[f][j][k][2];

                            double flux[5];

                            flux[0] = rb*(vb1*n1 + vb2*n2 + vb3*n3);
                            flux[1] = rb*vb1*(vb1*n1 + vb2*n2 + vb3*n3)+pb*n1;
                            flux[2] = rb*vb2*(vb1*n1 + vb2*n2 + vb3*n3)+pb*n2;
                            flux[3] = rb*vb3*(vb1*n1 + vb2*n2 + vb3*n3)+pb*n3;
                            flux[4] = (rEb+pb)*(vb1*n1 + vb2*n2 + vb3*n3);

                            mesh[block_id].xfluxes[f][j][k][0] = flux[0];
                            mesh[block_id].xfluxes[f][j][k][1] = flux[1];
                            mesh[block_id].xfluxes[f][j][k][2] = flux[2];
                            mesh[block_id].xfluxes[f][j][k][3] = flux[3];
                            mesh[block_id].xfluxes[f][j][k][4] = flux[4];

                        }
                    
                    }
                
                }
            
                if (bc==1){
                
                    for(int j=0; j<mesh[block_id].ny-1; j++){
                    
                        for(int k=0; k<mesh[block_id].nz-1; k++){
                        
                            int c = mesh[block_id].nx-1;
                            int f = c+1;

                            double rb  = mesh[block_id].states[c][j][k][0];
                            double vb1 = mesh[block_id].states[c][j][k][1]/rb;
                            double vb2 = mesh[block_id].states[c][j][k][2]/rb;
                            double vb3 = mesh[block_id].states[c][j][k][3]/rb;
                            double rEb = mesh[block_id].states[c][j][k][4];
                            double pb  = 0.4*(rEb - 0.5*rb*(vb1*vb1+vb2*vb2+vb3*vb3));
                            
                            double n1 = mesh[block_id].xn[f][j][k][0];
                            double n2 = mesh[block_id].xn[f][j][k][1];
                            double n3 = mesh[block_id].xn[f][j][k][2];

                            double flux[5];

                            flux[0] = rb*(vb1*n1 + vb2*n2 + vb3*n3);
                            flux[1] = rb*vb1*(vb1*n1 + vb2*n2 + vb3*n3)+pb*n1;
                            flux[2] = rb*vb2*(vb1*n1 + vb2*n2 + vb3*n3)+pb*n2;
                            flux[3] = rb*vb3*(vb1*n1 + vb2*n2 + vb3*n3)+pb*n3;
                            flux[4] = (rEb+pb)*(vb1*n1 + vb2*n2 + vb3*n3);

                            mesh[block_id].xfluxes[f][j][k][0] = flux[0];
                            mesh[block_id].xfluxes[f][j][k][1] = flux[1];
                            mesh[block_id].xfluxes[f][j][k][2] = flux[2];
                            mesh[block_id].xfluxes[f][j][k][3] = flux[3];
                            mesh[block_id].xfluxes[f][j][k][4] = flux[4];

                        }
                    
                    }
                
                }
            
                if (bc==2){
                
                    for(int i=0; i<mesh[block_id].nx-1; i++){
                    
                        for(int k=0; k<mesh[block_id].nz-1; k++){
                        
                            int c = 0;
                            int f = 0;

                            double rb  = mesh[block_id].states[i][c][k][0];
                            double vb1 = mesh[block_id].states[i][c][k][1]/rb;
                            double vb2 = mesh[block_id].states[i][c][k][2]/rb;
                            double vb3 = mesh[block_id].states[i][c][k][3]/rb;
                            double rEb = mesh[block_id].states[i][c][k][4];
                            double pb  = 0.4*(rEb - 0.5*rb*(vb1*vb1+vb2*vb2+vb3*vb3));
                            
                            double n1 = mesh[block_id].yn[i][f][k][0];
                            double n2 = mesh[block_id].yn[i][f][k][1];
                            double n3 = mesh[block_id].yn[i][f][k][2];

                            double flux[5];

                            flux[0] = rb*(vb1*n1 + vb2*n2 + vb3*n3);
                            flux[1] = rb*vb1*(vb1*n1 + vb2*n2 + vb3*n3)+pb*n1;
                            flux[2] = rb*vb2*(vb1*n1 + vb2*n2 + vb3*n3)+pb*n2;
                            flux[3] = rb*vb3*(vb1*n1 + vb2*n2 + vb3*n3)+pb*n3;
                            flux[4] = (rEb+pb)*(vb1*n1 + vb2*n2 + vb3*n3);

                            mesh[block_id].yfluxes[i][f][k][0] = flux[0];
                            mesh[block_id].yfluxes[i][f][k][1] = flux[1];
                            mesh[block_id].yfluxes[i][f][k][2] = flux[2];
                            mesh[block_id].yfluxes[i][f][k][3] = flux[3];
                            mesh[block_id].yfluxes[i][f][k][4] = flux[4];

                        }
                    
                    }
                
                }
            
                if (bc==3){
                
                    for(int i=0; i<mesh[block_id].nx-1; i++){
                    
                        for(int k=0; k<mesh[block_id].nz-1; k++){
                        
                            int c = mesh[block_id].ny-1;
                            int f = c+1;

                            double rb  = mesh[block_id].states[i][c][k][0];
                            double vb1 = mesh[block_id].states[i][c][k][1]/rb;
                            double vb2 = mesh[block_id].states[i][c][k][2]/rb;
                            double vb3 = mesh[block_id].states[i][c][k][3]/rb;
                            double rEb = mesh[block_id].states[i][c][k][4];
                            double pb  = 0.4*(rEb - 0.5*rb*(vb1*vb1+vb2*vb2+vb3*vb3));
                            
                            double n1 = mesh[block_id].yn[i][f][k][0];
                            double n2 = mesh[block_id].yn[i][f][k][1];
                            double n3 = mesh[block_id].yn[i][f][k][2];

                            double flux[5];

                            flux[0] = rb*(vb1*n1 + vb2*n2 + vb3*n3);
                            flux[1] = rb*vb1*(vb1*n1 + vb2*n2 + vb3*n3)+pb*n1;
                            flux[2] = rb*vb2*(vb1*n1 + vb2*n2 + vb3*n3)+pb*n2;
                            flux[3] = rb*vb3*(vb1*n1 + vb2*n2 + vb3*n3)+pb*n3;
                            flux[4] = (rEb+pb)*(vb1*n1 + vb2*n2 + vb3*n3);

                            mesh[block_id].yfluxes[i][f][k][0] = flux[0];
                            mesh[block_id].yfluxes[i][f][k][1] = flux[1];
                            mesh[block_id].yfluxes[i][f][k][2] = flux[2];
                            mesh[block_id].yfluxes[i][f][k][3] = flux[3];
                            mesh[block_id].yfluxes[i][f][k][4] = flux[4];

                        }
                    
                    }
                
                }
            
                if (bc==4){
                
                    for(int i=0; i<mesh[block_id].nx-1; i++){
                    
                        for(int j=0; j<mesh[block_id].ny-1; j++){
                        
                            int c = 0;
                            int f = 0;

                            double rb  = mesh[block_id].states[i][j][c][0];
                            double vb1 = mesh[block_id].states[i][j][c][1]/rb;
                            double vb2 = mesh[block_id].states[i][j][c][2]/rb;
                            double vb3 = mesh[block_id].states[i][j][c][3]/rb;
                            double rEb = mesh[block_id].states[i][j][c][4];
                            double pb  = 0.4*(rEb - 0.5*rb*(vb1*vb1+vb2*vb2+vb3*vb3));
                            
                            double n1 = mesh[block_id].zn[i][j][f][0];
                            double n2 = mesh[block_id].zn[i][j][f][1];
                            double n3 = mesh[block_id].zn[i][j][f][2];

                            double flux[5];

                            flux[0] = rb*(vb1*n1 + vb2*n2 + vb3*n3);
                            flux[1] = rb*vb1*(vb1*n1 + vb2*n2 + vb3*n3)+pb*n1;
                            flux[2] = rb*vb2*(vb1*n1 + vb2*n2 + vb3*n3)+pb*n2;
                            flux[3] = rb*vb3*(vb1*n1 + vb2*n2 + vb3*n3)+pb*n3;
                            flux[4] = (rEb+pb)*(vb1*n1 + vb2*n2 + vb3*n3);

                            mesh[block_id].zfluxes[i][j][f][0] = flux[0];
                            mesh[block_id].zfluxes[i][j][f][1] = flux[1];
                            mesh[block_id].zfluxes[i][j][f][2] = flux[2];
                            mesh[block_id].zfluxes[i][j][f][3] = flux[3];
                            mesh[block_id].zfluxes[i][j][f][4] = flux[4];

                        }
                    
                    }
                
                }
            
                if (bc==5){
                
                    for(int i=0; i<mesh[block_id].nx-1; i++){
                    
                        for(int j=0; j<mesh[block_id].ny-1; j++){
                        
                            int c = mesh[block_id].nz-1;
                            int f = c+1;

                            double rb  = mesh[block_id].states[i][j][c][0];
                            double vb1 = mesh[block_id].states[i][j][c][1]/rb;
                            double vb2 = mesh[block_id].states[i][j][c][2]/rb;
                            double vb3 = mesh[block_id].states[i][j][c][3]/rb;
                            double rEb = mesh[block_id].states[i][j][c][4];
                            double pb  = 0.4*(rEb - 0.5*rb*(vb1*vb1+vb2*vb2+vb3*vb3));
                            
                            double n1 = mesh[block_id].zn[i][j][f][0];
                            double n2 = mesh[block_id].zn[i][j][f][1];
                            double n3 = mesh[block_id].zn[i][j][f][2];

                            double flux[5];

                            flux[0] = rb*(vb1*n1 + vb2*n2 + vb3*n3);
                            flux[1] = rb*vb1*(vb1*n1 + vb2*n2 + vb3*n3)+pb*n1;
                            flux[2] = rb*vb2*(vb1*n1 + vb2*n2 + vb3*n3)+pb*n2;
                            flux[3] = rb*vb3*(vb1*n1 + vb2*n2 + vb3*n3)+pb*n3;
                            flux[4] = (rEb+pb)*(vb1*n1 + vb2*n2 + vb3*n3);

                            mesh[block_id].zfluxes[i][j][f][0] = flux[0];
                            mesh[block_id].zfluxes[i][j][f][1] = flux[1];
                            mesh[block_id].zfluxes[i][j][f][2] = flux[2];
                            mesh[block_id].zfluxes[i][j][f][3] = flux[3];
                            mesh[block_id].zfluxes[i][j][f][4] = flux[4];

                        }
                    
                    }
                
                }

            }
        
        }        
    
    }

}



void gather_fluxes(int nblock, block* mesh){

    bc_flux(nblock, mesh);

    for(int block_id=0; block_id<nblock; block_id++){
    

        for(int i=1; i<mesh[block_id].nx-1; i++){
        
            for(int j=0; j<mesh[block_id].ny-1; j++){
            
                for(int k=0; k<mesh[block_id].nz-1; k++){
                
                    Roe_flux(mesh[block_id].states[i-1][j][k],\
                             mesh[block_id].states[i][j][k],\
                             mesh[block_id].xn[i][j][k],\
                             mesh[block_id].xfluxes[i][j][k]);
                
                }
            
            }
        
        }        
    
        for(int i=0; i<mesh[block_id].nx-1; i++){
        
            for(int j=1; j<mesh[block_id].ny-1; j++){
            
                for(int k=0; k<mesh[block_id].nz-1; k++){
                
                    Roe_flux(mesh[block_id].states[i][j-1][k],\
                             mesh[block_id].states[i][j][k],\
                             mesh[block_id].yn[i][j][k],\
                             mesh[block_id].yfluxes[i][j][k]);
                
                }
            
            }
        
        }        
    
        for(int i=0; i<mesh[block_id].nx-1; i++){
        
            for(int j=0; j<mesh[block_id].ny-1; j++){
            
                for(int k=1; k<mesh[block_id].nz-1; k++){
                
                    Roe_flux(mesh[block_id].states[i][j][k-1],\
                             mesh[block_id].states[i][j][k],\
                             mesh[block_id].zn[i][j][k],\
                             mesh[block_id].zfluxes[i][j][k]);
                
                }
            
            }
        
        }

    
    }

}



void calc_residuals(int nblock, block* mesh){

    for(int block_id=0; block_id<nblock; block_id++){
    
        gather_fluxes(nblock,mesh);
        
        for(int i=0; i<mesh[block_id].nx-1; i++){
        
            for(int j=0; j<mesh[block_id].ny-1; j++){
            
                for(int k=0; k<mesh[block_id].nz-1; k++){
                
                    for(int l=0; l<5; l++){
                    
                        mesh[block_id].res[i][j][k][l] = mesh[block_id].xfluxes[i][j][k][l]+\
                                                         mesh[block_id].yfluxes[i][j][k][l]+\
                                                         mesh[block_id].zfluxes[i][j][k][l]-\
                                                         mesh[block_id].xfluxes[i+1][j][k][l]-\
                                                         mesh[block_id].yfluxes[i][j+1][k][l]-\
                                                         mesh[block_id].zfluxes[i][j][k+1][l];
                    
                    }
                
                }
            
            }
        
        }
    
    }

}
