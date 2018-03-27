#include <iostream>
#include <fstream>
#include <cmath>
#include "var_struct.cpp"

using namespace std;

int read_mesh(char* meshfile, block* mesh){

    int nblock;
    double s1x, s2x, s3x, s4x,
           s1y, s2y, s3y, s4y,
           s1z, s2z, s3z, s4z;

    ifstream mfile;
    mfile.open(meshfile);

    mfile>>nblock;
    mesh = new block[nblock];

    cout<<"============================================================================================="<<endl;
    cout<<"Setting up mesh, boundary conditions, states and fluxes"<<endl;
    cout<<"============================================================================================="<<endl;
    cout<<endl;

    for(int block_id=0; block_id<nblock; block_id++){

        cout<<"Initializing block "<<block_id+1<<endl;
    
        mfile>>mesh[block_id].nx>>mesh[block_id].ny>>mesh[block_id].nz;

        
        
        // Initialize coordinates

        mesh[block_id].x  =   new double**[mesh[block_id].nx];
        mesh[block_id].y  =   new double**[mesh[block_id].nx];
        mesh[block_id].z  =   new double**[mesh[block_id].nx];

        for(int i=0; i<mesh[block_id].nx; i++){
        
            mesh[block_id].x[i]  =   new double*[mesh[block_id].ny];
            mesh[block_id].y[i]  =   new double*[mesh[block_id].ny];
            mesh[block_id].z[i]  =   new double*[mesh[block_id].ny];

            for(int j=0; j<mesh[block_id].ny; j++){
            
                mesh[block_id].x[i][j] = new double[mesh[block_id].nz];
                mesh[block_id].y[i][j] = new double[mesh[block_id].nz];
                mesh[block_id].z[i][j] = new double[mesh[block_id].nz];
            
            }
        
        }



        // Initialize x-normals, x-faces and x-fluxes
    
        mesh[block_id].xn       =   new double***[mesh[block_id].nx];
        mesh[block_id].xfc      =   new double***[mesh[block_id].nx];
        mesh[block_id].xfluxes  =   new double***[mesh[block_id].nx];
        mesh[block_id].xArea    =   new double**[mesh[block_id].nx];

        for(int i=0; i<mesh[block_id].nx; i++){
        
            mesh[block_id].xn[i]        =   new double**[mesh[block_id].ny-1];
            mesh[block_id].xfc[i]       =   new double**[mesh[block_id].ny-1];
            mesh[block_id].xfluxes[i]   =   new double**[mesh[block_id].ny-1];
            mesh[block_id].xArea[i]     =   new double*[mesh[block_id].ny-1];

            for(int j=0; j<mesh[block_id].ny-1; j++){
            
                mesh[block_id].xn[i][j]         = new double*[mesh[block_id].nz-1];
                mesh[block_id].xfc[i][j]        = new double*[mesh[block_id].nz-1];
                mesh[block_id].xfluxes[i][j]    = new double*[mesh[block_id].nz-1];
                mesh[block_id].xArea[i][j]      = new double[mesh[block_id].nz-1];

                for(int k=0; k<mesh[block_id].nz-1; k++){
                
                    mesh[block_id].xn[i][j][k]      = new double[3];
                    mesh[block_id].xfc[i][j][k]     = new double[3];
                    mesh[block_id].xfluxes[i][j][k] = new double[5];
                
                }
            
            }
        
        }



        // Initialize y-normals, y-faces and y-fluxes

        mesh[block_id].yn       =   new double***[mesh[block_id].nx-1];
        mesh[block_id].yfc      =   new double***[mesh[block_id].nx-1];
        mesh[block_id].yfluxes  =   new double***[mesh[block_id].nx-1];
        mesh[block_id].yArea    =   new double**[mesh[block_id].nx-1];

        for(int i=0; i<mesh[block_id].nx-1; i++){
        
            mesh[block_id].yn[i]        =   new double**[mesh[block_id].ny];
            mesh[block_id].yfc[i]       =   new double**[mesh[block_id].ny];
            mesh[block_id].yfluxes[i]   =   new double**[mesh[block_id].ny];
            mesh[block_id].yArea[i]     =   new double*[mesh[block_id].ny];

            for(int j=0; j<mesh[block_id].ny; j++){
            
                mesh[block_id].yn[i][j]         = new double*[mesh[block_id].nz-1];
                mesh[block_id].yfc[i][j]        = new double*[mesh[block_id].nz-1];
                mesh[block_id].yfluxes[i][j]    = new double*[mesh[block_id].nz-1];
                mesh[block_id].yArea[i][j]      = new double[mesh[block_id].nz-1];
            
                for(int k=0; k<mesh[block_id].nz-1; k++){
                
                    mesh[block_id].yn[i][j][k]      = new double[3];
                    mesh[block_id].yfc[i][j][k]     = new double[3];
                    mesh[block_id].yfluxes[i][j][k] = new double[5];
                
                }

            }
        
        }



        // Initialize z-normals, z-faces and z-fluxes

        mesh[block_id].zn       =   new double***[mesh[block_id].nx-1];
        mesh[block_id].zfc      =   new double***[mesh[block_id].nx-1];
        mesh[block_id].zfluxes  =   new double***[mesh[block_id].nx-1];
        mesh[block_id].zArea    =   new double**[mesh[block_id].nx-1];

        for(int i=0; i<mesh[block_id].nx-1; i++){
        
            mesh[block_id].zn[i]        =   new double**[mesh[block_id].ny-1];
            mesh[block_id].zfc[i]       =   new double**[mesh[block_id].ny-1];
            mesh[block_id].zfluxes[i]   =   new double**[mesh[block_id].ny-1];
            mesh[block_id].zArea[i]     =   new double*[mesh[block_id].ny-1];

            for(int j=0; j<mesh[block_id].ny-1; j++){
            
                mesh[block_id].zn[i][j]         = new double*[mesh[block_id].nz];
                mesh[block_id].zfc[i][j]        = new double*[mesh[block_id].nz];
                mesh[block_id].zfluxes[i][j]    = new double*[mesh[block_id].nz];
                mesh[block_id].zArea[i][j]      = new double[mesh[block_id].nz];
            
                for(int k=0; k<mesh[block_id].nz; k++){
                
                    mesh[block_id].zn[i][j][k]      = new double[3];
                    mesh[block_id].zfc[i][j][k]     = new double[3];
                    mesh[block_id].zfluxes[i][j][k] = new double[5];
                
                }

            }
        
        }



        // Initialize cell-centroids and states
        
        mesh[block_id].cc       =   new double***[mesh[block_id].nx-1];
        mesh[block_id].states   =   new double***[mesh[block_id].nx-1];
        mesh[block_id].res      =   new double***[mesh[block_id].nx-1];
        mesh[block_id].tempres  =   new double***[mesh[block_id].nx-1];
        mesh[block_id].Volume   =   new double**[mesh[block_id].nx-1];

        for(int i=0; i<mesh[block_id].nx-1; i++){
        
            mesh[block_id].cc[i]        =   new double**[mesh[block_id].ny-1];
            mesh[block_id].states[i]    =   new double**[mesh[block_id].ny-1];
            mesh[block_id].res[i]       =   new double**[mesh[block_id].ny-1];
            mesh[block_id].tempres[i]   =   new double**[mesh[block_id].ny-1];
            mesh[block_id].Volume[i]    =   new double*[mesh[block_id].ny-1];

            for(int j=0; j<mesh[block_id].ny-1; j++){
            
                mesh[block_id].cc[i][j]     = new double*[mesh[block_id].nz-1];
                mesh[block_id].states[i][j] = new double*[mesh[block_id].nz-1];
                mesh[block_id].res[i][j]    = new double*[mesh[block_id].nz-1];
                mesh[block_id].tempres[i][j]= new double*[mesh[block_id].nz-1];
                mesh[block_id].Volume[i][j] = new double[mesh[block_id].nz-1];

                for(int k=0; k<mesh[block_id].nz-1; k++){
                
                    mesh[block_id].cc[i][j][k]      = new double[3];
                    mesh[block_id].states[i][j][k]  = new double[5];  
                    mesh[block_id].res[i][j][k]     = new double[5];  
                    mesh[block_id].tempres[i][j][k] = new double[5];  
                    mesh[block_id].Volume[i][j][k]  = 0.0;
                
                }
            
            }
        
        }
        


        cout<<"Gathering coordinates..."<<endl;

        // Assign coordinates

        for(int k=0; k<mesh[block_id].nz; k++){
        
            for(int j=0; j<mesh[block_id].ny; j++){
            
                for(int i=0; i<mesh[block_id].nx; i++){
                
                    mfile>>mesh[block_id].x[i][j][k]>>mesh[block_id].y[i][j][k]>>mesh[block_id].z[i][j][k];
                
                }
            
            }
        
        }



        cout<<"Setting up boundary conditions..."<<endl;

        // Assign block comms

        for(int i=0; i<6; i++){
        
            mfile>>mesh[block_id].bc_type[i];
        
        }



        // Assign boundary states
        
        for(int i=0; i<6; i++){

            for(int j=0; j<5; j++){
        
                mfile>>mesh[block_id].bc_states[i][j];
        
            }

        }



        cout<<"Calculating ξ-face geometry parameters..."<<endl;

        // Calculate face centroid, Area and normal for x-faces
        
        for(int i=0; i<mesh[block_id].nx; i++){
        
            for(int j=0; j<mesh[block_id].ny-1; j++){
            
                for(int k=0; k<mesh[block_id].nz-1; k++){
                
                    mesh[block_id].xfc[i][j][k][0] = 0.25*(mesh[block_id].x[i][j][k]+\
                                                           mesh[block_id].x[i][j+1][k]+\
                                                           mesh[block_id].x[i][j][k+1]+\
                                                           mesh[block_id].x[i][j+1][k+1]);
                
                    mesh[block_id].xfc[i][j][k][1] = 0.25*(mesh[block_id].y[i][j][k]+\
                                                           mesh[block_id].y[i][j+1][k]+\
                                                           mesh[block_id].y[i][j][k+1]+\
                                                           mesh[block_id].y[i][j+1][k+1]);
                
                    mesh[block_id].xfc[i][j][k][2] = 0.25*(mesh[block_id].z[i][j][k]+\
                                                           mesh[block_id].z[i][j+1][k]+\
                                                           mesh[block_id].z[i][j][k+1]+\
                                                           mesh[block_id].z[i][j+1][k+1]);

                    s1x   =   mesh[block_id].x[i][j+1][k]  -mesh[block_id].x[i][j][k];
                    s2x   =   mesh[block_id].x[i][j+1][k+1]-mesh[block_id].x[i][j+1][k];
                    s3x   =   mesh[block_id].x[i][j][k+1]  -mesh[block_id].x[i][j+1][k+1];
                    s4x   =   mesh[block_id].x[i][j][k]    -mesh[block_id].x[i][j][k+1];

                    s1y   =   mesh[block_id].y[i][j+1][k]  -mesh[block_id].y[i][j][k];
                    s2y   =   mesh[block_id].y[i][j+1][k+1]-mesh[block_id].y[i][j+1][k];
                    s3y   =   mesh[block_id].y[i][j][k+1]  -mesh[block_id].y[i][j+1][k+1];
                    s4y   =   mesh[block_id].y[i][j][k]    -mesh[block_id].y[i][j][k+1];

                    s1z   =   mesh[block_id].z[i][j+1][k]  -mesh[block_id].z[i][j][k];
                    s2z   =   mesh[block_id].z[i][j+1][k+1]-mesh[block_id].z[i][j+1][k];
                    s3z   =   mesh[block_id].z[i][j][k+1]  -mesh[block_id].z[i][j+1][k+1];
                    s4z   =   mesh[block_id].z[i][j][k]    -mesh[block_id].z[i][j][k+1];

                    mesh[block_id].xn[i][j][k][0]  = 0.5*(s1y*s2z-s1z*s2y);
                    mesh[block_id].xn[i][j][k][1]  = 0.5*(s1z*s2x-s1x*s2z);
                    mesh[block_id].xn[i][j][k][2]  = 0.5*(s1x*s2y-s1y*s2x);
               
                    mesh[block_id].xArea[i][j][k]  = sqrt(pow(mesh[block_id].xn[i][j][k][0],2)+\
                                                          pow(mesh[block_id].xn[i][j][k][1],2)+\
                                                          pow(mesh[block_id].xn[i][j][k][2],2));
                    
                    mesh[block_id].xn[i][j][k][0]  = mesh[block_id].xn[i][j][k][0]/\
                                                     mesh[block_id].xArea[i][j][k];
                    mesh[block_id].xn[i][j][k][1]  = mesh[block_id].xn[i][j][k][1]/\
                                                     mesh[block_id].xArea[i][j][k];
                    mesh[block_id].xn[i][j][k][2]  = mesh[block_id].xn[i][j][k][2]/\
                                                     mesh[block_id].xArea[i][j][k];

                }
            
            }
        
        }



        cout<<"Calculating η-face geometry parameters..."<<endl;

        // Calculate face centroid, Area and normal for y-faces
        
        for(int i=0; i<mesh[block_id].nx-1; i++){
        
            for(int j=0; j<mesh[block_id].ny; j++){
            
                for(int k=0; k<mesh[block_id].nz-1; k++){
                
                    mesh[block_id].yfc[i][j][k][0] = 0.25*(mesh[block_id].x[i][j][k]+\
                                                           mesh[block_id].x[i+1][j][k]+\
                                                           mesh[block_id].x[i][j][k+1]+\
                                                           mesh[block_id].x[i+1][j][k+1]);
                
                    mesh[block_id].yfc[i][j][k][1] = 0.25*(mesh[block_id].y[i][j][k]+\
                                                           mesh[block_id].y[i+1][j][k]+\
                                                           mesh[block_id].y[i][j][k+1]+\
                                                           mesh[block_id].y[i+1][j][k+1]);
                
                    mesh[block_id].yfc[i][j][k][2] = 0.25*(mesh[block_id].z[i][j][k]+\
                                                           mesh[block_id].z[i+1][j][k]+\
                                                           mesh[block_id].z[i][j][k+1]+\
                                                           mesh[block_id].z[i+1][j][k+1]);

                    s1x   =   mesh[block_id].x[i][j][k+1]  -mesh[block_id].x[i][j][k];
                    s2x   =   mesh[block_id].x[i+1][j][k+1]-mesh[block_id].x[i][j][k+1];
                    s3x   =   mesh[block_id].x[i+1][j][k]  -mesh[block_id].x[i+1][j][k+1];
                    s4x   =   mesh[block_id].x[i][j][k]    -mesh[block_id].x[i+1][j][k];

                    s1y   =   mesh[block_id].y[i][j][k+1]  -mesh[block_id].y[i][j][k];
                    s2y   =   mesh[block_id].y[i+1][j][k+1]-mesh[block_id].y[i][j][k+1];
                    s3y   =   mesh[block_id].y[i+1][j][k]  -mesh[block_id].y[i+1][j][k+1];
                    s4y   =   mesh[block_id].y[i][j][k]    -mesh[block_id].y[i+1][j][k];

                    s1z   =   mesh[block_id].z[i][j][k+1]  -mesh[block_id].z[i][j][k];
                    s2z   =   mesh[block_id].z[i+1][j][k+1]-mesh[block_id].z[i][j][k+1];
                    s3z   =   mesh[block_id].z[i+1][j][k]  -mesh[block_id].z[i+1][j][k+1];
                    s4z   =   mesh[block_id].z[i][j][k]    -mesh[block_id].z[i+1][j][k];

                    mesh[block_id].yn[i][j][k][0]  = 0.5*(s1y*s2z-s1z*s2y);
                    mesh[block_id].yn[i][j][k][1]  = 0.5*(s1z*s2x-s1x*s2z);
                    mesh[block_id].yn[i][j][k][2]  = 0.5*(s1x*s2y-s1y*s2x);
               
                    mesh[block_id].yArea[i][j][k]  = sqrt(pow(mesh[block_id].yn[i][j][k][0],2)+\
                                                          pow(mesh[block_id].yn[i][j][k][1],2)+\
                                                          pow(mesh[block_id].yn[i][j][k][2],2));
                    
                    mesh[block_id].yn[i][j][k][0]  = mesh[block_id].yn[i][j][k][0]/\
                                                     mesh[block_id].yArea[i][j][k];
                    mesh[block_id].yn[i][j][k][1]  = mesh[block_id].yn[i][j][k][1]/\
                                                     mesh[block_id].yArea[i][j][k];
                    mesh[block_id].yn[i][j][k][2]  = mesh[block_id].yn[i][j][k][2]/\
                                                     mesh[block_id].yArea[i][j][k];

                }
            
            }
        
        }



        cout<<"Calculating ζ-face geometry parameters..."<<endl;

        // Calculate face centroid, Area and normal for z-faces
        
        for(int i=0; i<mesh[block_id].nx-1; i++){
        
            for(int j=0; j<mesh[block_id].ny-1; j++){
            
                for(int k=0; k<mesh[block_id].nz; k++){
                
                    mesh[block_id].zfc[i][j][k][0] = 0.25*(mesh[block_id].x[i][j][k]+\
                                                           mesh[block_id].x[i][j+1][k]+\
                                                           mesh[block_id].x[i+1][j][k]+\
                                                           mesh[block_id].x[i+1][j+1][k]);
                
                    mesh[block_id].zfc[i][j][k][1] = 0.25*(mesh[block_id].y[i][j][k]+\
                                                           mesh[block_id].y[i][j+1][k]+\
                                                           mesh[block_id].y[i+1][j][k]+\
                                                           mesh[block_id].y[i+1][j+1][k]);
                
                    mesh[block_id].zfc[i][j][k][2] = 0.25*(mesh[block_id].z[i][j][k]+\
                                                           mesh[block_id].z[i][j+1][k]+\
                                                           mesh[block_id].z[i+1][j][k]+\
                                                           mesh[block_id].z[i+1][j+1][k]);

                    s1x   =   mesh[block_id].x[i+1][j][k]  -mesh[block_id].x[i][j][k];
                    s2x   =   mesh[block_id].x[i+1][j+1][k]-mesh[block_id].x[i+1][j][k];
                    s3x   =   mesh[block_id].x[i][j+1][k]  -mesh[block_id].x[i+1][j+1][k];
                    s4x   =   mesh[block_id].x[i][j][k]    -mesh[block_id].x[i][j+1][k];

                    s1y   =   mesh[block_id].y[i+1][j][k]  -mesh[block_id].y[i][j][k];
                    s2y   =   mesh[block_id].y[i+1][j+1][k]-mesh[block_id].y[i+1][j][k];
                    s3y   =   mesh[block_id].y[i][j+1][k]  -mesh[block_id].y[i+1][j+1][k];
                    s4y   =   mesh[block_id].y[i][j][k]    -mesh[block_id].y[i][j+1][k];

                    s1z   =   mesh[block_id].z[i+1][j][k]  -mesh[block_id].z[i][j][k];
                    s2z   =   mesh[block_id].z[i+1][j+1][k]-mesh[block_id].z[i+1][j][k];
                    s3z   =   mesh[block_id].z[i][j+1][k]  -mesh[block_id].z[i+1][j+1][k];
                    s4z   =   mesh[block_id].z[i][j][k]    -mesh[block_id].z[i][j+1][k];

                    mesh[block_id].zn[i][j][k][0]  = 0.5*(s1y*s2z-s1z*s2y);
                    mesh[block_id].zn[i][j][k][1]  = 0.5*(s1z*s2x-s1x*s2z);
                    mesh[block_id].zn[i][j][k][2]  = 0.5*(s1x*s2y-s1y*s2x);
               
                    mesh[block_id].zArea[i][j][k]  = sqrt(pow(mesh[block_id].zn[i][j][k][0],2)+\
                                                          pow(mesh[block_id].zn[i][j][k][1],2)+\
                                                          pow(mesh[block_id].zn[i][j][k][2],2));
                    
                    mesh[block_id].zn[i][j][k][0]  = mesh[block_id].zn[i][j][k][0]/\
                                                     mesh[block_id].zArea[i][j][k];
                    mesh[block_id].zn[i][j][k][1]  = mesh[block_id].zn[i][j][k][1]/\
                                                     mesh[block_id].zArea[i][j][k];
                    mesh[block_id].zn[i][j][k][2]  = mesh[block_id].zn[i][j][k][2]/\
                                                     mesh[block_id].zArea[i][j][k];

                }
            
            }
        
        }

        
        
        cout<<"Calculating cell geometry parameters..."<<endl;

        // Calculate cell centroid and cell volume
        
        for(int i=0; i<mesh[block_id].nx-1; i++){
        
            for(int j=0; j<mesh[block_id].ny-1; j++){
            
                for(int k=0; k<mesh[block_id].nz-1; k++){
                
                    mesh[block_id].cc[i][j][k][0]  = 0.125*(mesh[block_id].x[i][j][k]+\
                                                            mesh[block_id].x[i+1][j][k]+\
                                                            mesh[block_id].x[i][j+1][k]+\
                                                            mesh[block_id].x[i][j][k+1]+\
                                                            mesh[block_id].x[i+1][j+1][k]+\
                                                            mesh[block_id].x[i+1][j][k+1]+\
                                                            mesh[block_id].x[i][j+1][k+1]+\
                                                            mesh[block_id].x[i+1][j+1][k+1]);
                
                    mesh[block_id].cc[i][j][k][1]  = 0.125*(mesh[block_id].y[i][j][k]+\
                                                            mesh[block_id].y[i+1][j][k]+\
                                                            mesh[block_id].y[i][j+1][k]+\
                                                            mesh[block_id].y[i][j][k+1]+\
                                                            mesh[block_id].y[i+1][j+1][k]+\
                                                            mesh[block_id].y[i+1][j][k+1]+\
                                                            mesh[block_id].y[i][j+1][k+1]+\
                                                            mesh[block_id].y[i+1][j+1][k+1]);
                
                    mesh[block_id].cc[i][j][k][2]  = 0.125*(mesh[block_id].z[i][j][k]+\
                                                            mesh[block_id].z[i+1][j][k]+\
                                                            mesh[block_id].z[i][j+1][k]+\
                                                            mesh[block_id].z[i][j][k+1]+\
                                                            mesh[block_id].z[i+1][j+1][k]+\
                                                            mesh[block_id].z[i+1][j][k+1]+\
                                                            mesh[block_id].z[i][j+1][k+1]+\
                                                            mesh[block_id].z[i+1][j+1][k+1]);
                
                    
                    
                    s1x   =   mesh[block_id].x[i+1][j][k] - mesh[block_id].x[i][j][k];
                    s2x   =   mesh[block_id].x[i][j+1][k] - mesh[block_id].x[i][j][k];
                    s3x   =   mesh[block_id].x[i][j][k+1] - mesh[block_id].x[i][j][k];

                    s1y   =   mesh[block_id].y[i+1][j][k] - mesh[block_id].y[i][j][k];
                    s2y   =   mesh[block_id].y[i][j+1][k] - mesh[block_id].y[i][j][k];
                    s3y   =   mesh[block_id].y[i][j][k+1] - mesh[block_id].y[i][j][k];

                    s1z   =   mesh[block_id].z[i+1][j][k] - mesh[block_id].z[i][j][k];
                    s2z   =   mesh[block_id].z[i][j+1][k] - mesh[block_id].z[i][j][k];
                    s3z   =   mesh[block_id].z[i][j][k+1] - mesh[block_id].z[i][j][k];

                    mesh[block_id].Volume[i][j][k] = mesh[block_id].Volume[i][j][k]+\
                                                     1/6*abs(s1x*(s2y*s3z-s2z*s3y)+\
                                                             s1y*(s2z*s3x-s2x*s3z)+\
                                                             s1z*(s2x*s3y-s2y*s3x));

                    
                    
                    s1x   =   mesh[block_id].x[i+1][j][k] - mesh[block_id].x[i+1][j+1][k];
                    s2x   =   mesh[block_id].x[i][j+1][k] - mesh[block_id].x[i+1][j+1][k];
                    s3x   =   mesh[block_id].x[i+1][j+1][k+1] - mesh[block_id].x[i+1][j+1][k];

                    s1y   =   mesh[block_id].y[i+1][j][k] - mesh[block_id].y[i+1][j+1][k];
                    s2y   =   mesh[block_id].y[i][j+1][k] - mesh[block_id].y[i+1][j+1][k];
                    s3y   =   mesh[block_id].y[i+1][j+1][k+1] - mesh[block_id].y[i+1][j+1][k];

                    s1z   =   mesh[block_id].z[i+1][j][k] - mesh[block_id].z[i+1][j+1][k];
                    s2z   =   mesh[block_id].z[i][j+1][k] - mesh[block_id].z[i+1][j+1][k];
                    s3z   =   mesh[block_id].z[i+1][j+1][k+1] - mesh[block_id].z[i+1][j+1][k];

                    mesh[block_id].Volume[i][j][k] = mesh[block_id].Volume[i][j][k]+\
                                                     1/6*abs(s1x*(s2y*s3z-s2z*s3y)+\
                                                             s1y*(s2z*s3x-s2x*s3z)+\
                                                             s1z*(s2x*s3y-s2y*s3x));

                    
                    
                    s1x   =   mesh[block_id].x[i+1][j+1][k+1] - mesh[block_id].x[i][j+1][k+1];
                    s2x   =   mesh[block_id].x[i][j+1][k] - mesh[block_id].x[i][j+1][k+1];
                    s3x   =   mesh[block_id].x[i][j][k+1] - mesh[block_id].x[i][j+1][k+1];

                    s1y   =   mesh[block_id].y[i+1][j+1][k+1] - mesh[block_id].y[i][j+1][k+1];
                    s2y   =   mesh[block_id].y[i][j+1][k] - mesh[block_id].y[i][j+1][k+1];
                    s3y   =   mesh[block_id].y[i][j][k+1] - mesh[block_id].y[i][j+1][k+1];

                    s1z   =   mesh[block_id].z[i+1][j+1][k+1] - mesh[block_id].z[i][j+1][k+1];
                    s2z   =   mesh[block_id].z[i][j+1][k] - mesh[block_id].z[i][j+1][k+1];
                    s3z   =   mesh[block_id].z[i][j][k+1] - mesh[block_id].z[i][j+1][k+1];

                    mesh[block_id].Volume[i][j][k] = mesh[block_id].Volume[i][j][k]+\
                                                     1/6*abs(s1x*(s2y*s3z-s2z*s3y)+\
                                                             s1y*(s2z*s3x-s2x*s3z)+\
                                                             s1z*(s2x*s3y-s2y*s3x));

                    
                    
                    s1x   =   mesh[block_id].x[i+1][j][k] - mesh[block_id].x[i+1][j][k+1];
                    s2x   =   mesh[block_id].x[i+1][j+1][k+1] - mesh[block_id].x[i+1][j][k+1];
                    s3x   =   mesh[block_id].x[i][j][k+1] - mesh[block_id].x[i+1][j][k+1];

                    s1y   =   mesh[block_id].y[i+1][j][k] - mesh[block_id].y[i+1][j][k+1];
                    s2y   =   mesh[block_id].y[i+1][j+1][k+1] - mesh[block_id].y[i+1][j][k+1];
                    s3y   =   mesh[block_id].y[i][j][k+1] - mesh[block_id].y[i+1][j][k+1];

                    s1z   =   mesh[block_id].z[i+1][j][k] - mesh[block_id].z[i+1][j][k+1];
                    s2z   =   mesh[block_id].z[i+1][j+1][k+1] - mesh[block_id].z[i+1][j][k+1];
                    s3z   =   mesh[block_id].z[i][j][k+1] - mesh[block_id].z[i+1][j][k+1];

                    mesh[block_id].Volume[i][j][k] = mesh[block_id].Volume[i][j][k]+\
                                                     1/6*abs(s1x*(s2y*s3z-s2z*s3y)+\
                                                             s1y*(s2z*s3x-s2x*s3z)+\
                                                             s1z*(s2x*s3y-s2y*s3x));

                    
                    
                    s1x   =   mesh[block_id].x[i+1][j][k] - mesh[block_id].x[i][j][k+1];
                    s2x   =   mesh[block_id].x[i][j+1][k] - mesh[block_id].x[i][j][k+1];
                    s3x   =   mesh[block_id].x[i+1][j+1][k+1] - mesh[block_id].x[i][j][k+1];

                    s1y   =   mesh[block_id].y[i+1][j][k] - mesh[block_id].y[i][j][k+1];
                    s2y   =   mesh[block_id].y[i][j+1][k] - mesh[block_id].y[i][j][k+1];
                    s3y   =   mesh[block_id].y[i+1][j+1][k+1] - mesh[block_id].y[i][j][k+1];

                    s1z   =   mesh[block_id].z[i+1][j][k] - mesh[block_id].z[i][j][k+1];
                    s2z   =   mesh[block_id].z[i][j+1][k] - mesh[block_id].z[i][j][k+1];
                    s3z   =   mesh[block_id].z[i+1][j+1][k+1] - mesh[block_id].z[i][j][k+1];

                    mesh[block_id].Volume[i][j][k] = mesh[block_id].Volume[i][j][k]+\
                                                     1/6*abs(s1x*(s2y*s3z-s2z*s3y)+\
                                                             s1y*(s2z*s3x-s2x*s3z)+\
                                                             s1z*(s2x*s3y-s2y*s3x));

                }
            
            }
        
        }

    }

    mfile.close();
    return nblock;

}
