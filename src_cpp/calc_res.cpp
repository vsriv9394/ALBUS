#include <iostream>
#include <cmath>
#include "flux.cpp"

using namespace std;



void bc_flux(int nblock, block* mesh, partition* map){

    for(int block_id=0; block_id<nblock; block_id++){
        
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
        
        // Sending and receiving data from left to the right cells
        if (mesh[block_id].bc_type[1][0]>=0){
        
            int proc_r      = mesh[block_id].bc_type[1][0];
            int proc_r_id   = mesh[block_id].bc_type[1][1];
            
            if(proc_r!=rank){

                int siz_s[4]  = {mesh[block_id].nx-1, mesh[block_id].ny-1, mesh[block_id].nz-1, 5};
                int siz_sa[4] = {1, mesh[block_id].ny-1, mesh[block_id].nz-1, 5};
                int starts[4]    = {mesh[block_id].nx-2,0,0,0};
                MPI_Datatype Right_boundary;
                MPI_Type_create_subarray(4, siz_s, siz_sa, starts, MPI_ORDER_C, MPI_DOUBLE, &Right_boundary);
                MPI_Type_commit(&Right_boundary);
                MPI_Request* req;
                MPI_Isend(&(mesh[block_id].states[0][0][0][0]), 1, Right_boundary, proc_r, proc_r*1000+proc_r_id, MPI_COMM_WORLD, req);
            
            }

        }

        if (mesh[block_id].bc_type[0][0]>=0){

            int proc_l      = mesh[block_id].bc_type[0][0];
            int proc_l_id   = mesh[block_id].bc_type[0][1];
            
            double**** left_states;
            left_states     = new double***[1];
            left_states[0]  = new double**[mesh[block_id].ny-1];
            
            for(int j=0; j<mesh[block_id].ny-1; j++){

                left_states[0][j] = new double*[mesh[block_id].nz-1];
                    
                for(int k=0; k<mesh[block_id].nz-1; k++){
            
                    left_states[0][j][k] = new double[5];
            
                }

            }

            if(proc_l!=rank){
            
                MPI_Request* req;
                MPI_Irecv(&(left_states[0][0][0][0]), (mesh[block_id].ny-1)*(mesh[block_id].nz-1)*5, MPI_DOUBLE, proc_l, rank*1000+block_id, MPI_COMM_WORLD, req);
                MPI_Wait(req, MPI_STATUS_IGNORE);
                
            }

            for(int j=0; j<mesh[block_id].ny-1; j++){
                    
                for(int k=0; k<mesh[block_id].nz-1; k++){
            
                    if(proc_l==rank){

                        for(int sid=0; sid<5; sid++){

                            left_states[0][j][k][sid] = mesh[proc_l_id].states[mesh[proc_l_id].nx-2][j][k][sid];

                        }

                    }

                    Roe_flux( left_states[0][j][k],\
                              mesh[ block_id ].states[0][j][k],\
                              mesh[ block_id ].xn[0][j][k],\
                              mesh[ block_id ].xfluxes[0][j][k]);
            
                }

            }

            delete[] left_states;
            
        }
        
        // Sending and receiving data from right to the left cells
        if (mesh[block_id].bc_type[0][0]>=0){
        
            int proc_l      = mesh[block_id].bc_type[0][0];
            int proc_l_id   = mesh[block_id].bc_type[0][1];
            
            if(proc_l!=rank){

                int siz_s[4]  = {mesh[block_id].nx-1, mesh[block_id].ny-1, mesh[block_id].nz-1, 5};
                int siz_sa[4] = {1, mesh[block_id].ny-1, mesh[block_id].nz-1, 5};
                int starts[4] = {0,0,0,0};
                MPI_Datatype Left_boundary;
                MPI_Type_create_subarray(4, siz_s, siz_sa, starts, MPI_ORDER_C, MPI_DOUBLE, &Left_boundary);
                MPI_Type_commit(&Left_boundary);
                MPI_Request* req;
                MPI_Isend(&(mesh[block_id].states[0][0][0][0]), 1, Left_boundary, proc_l, proc_l*1000+proc_l_id, MPI_COMM_WORLD, req);
            
            }

        }

        if (mesh[block_id].bc_type[1][0]>=0){

            int proc_r      = mesh[block_id].bc_type[1][0];
            int proc_r_id   = mesh[block_id].bc_type[1][1];
            
            double**** right_states;
            right_states     = new double***[1];
            right_states[0]  = new double**[mesh[block_id].ny-1];
            
            for(int j=0; j<mesh[block_id].ny-1; j++){

                right_states[0][j] = new double*[mesh[block_id].nz-1];
                    
                for(int k=0; k<mesh[block_id].nz-1; k++){
            
                    right_states[0][j][k] = new double[5];
            
                }

            }

            if(proc_r!=rank){
            
                MPI_Request* req;
                MPI_Irecv(&(right_states[0][0][0][0]), (mesh[block_id].ny-1)*(mesh[block_id].nz-1)*5, MPI_DOUBLE, proc_r, rank*1000+block_id, MPI_COMM_WORLD, req);
                MPI_Wait(req, MPI_STATUS_IGNORE);
                
            }

            for(int j=0; j<mesh[block_id].ny-1; j++){
                    
                for(int k=0; k<mesh[block_id].nz-1; k++){
            
                    if(proc_r==rank){

                        for(int sid=0; sid<5; sid++){

                            right_states[0][j][k][sid] = mesh[proc_r_id].states[0][j][k][sid];

                        }

                    }

                    Roe_flux( mesh[ block_id ].states[mesh[block_id].nx-2][j][k],\
                              right_states[0][j][k],\
                              mesh[ block_id ].xn[mesh[block_id].nx-1][j][k],\
                              mesh[ block_id ].xfluxes[mesh[block_id].nx-1][j][k]);
            
                }

            }

            delete[] right_states;
            
        }

        // Sending and receiving data from bottom to the top cells
        if (mesh[block_id].bc_type[3][0]>=0){
        
            int proc_r      = mesh[block_id].bc_type[3][0];
            int proc_r_id   = mesh[block_id].bc_type[3][1];
            
            if(proc_r!=rank){

                int siz_s[4]  = {mesh[block_id].nx-1, mesh[block_id].ny-1, mesh[block_id].nz-1, 5};
                int siz_sa[4] = {mesh[block_id].nx-1, 1, mesh[block_id].nz-1, 5};
                int starts[4]    = {0,mesh[block_id].ny-2,0,0};
                MPI_Datatype Right_boundary;
                MPI_Type_create_subarray(4, siz_s, siz_sa, starts, MPI_ORDER_C, MPI_DOUBLE, &Right_boundary);
                MPI_Type_commit(&Right_boundary);
                MPI_Request* req;
                MPI_Isend(&(mesh[block_id].states[0][0][0][0]), 1, Right_boundary, proc_r, proc_r*1000+proc_r_id, MPI_COMM_WORLD, req);
            
            }

        }

        if (mesh[block_id].bc_type[2][0]>=0){

            int proc_l      = mesh[block_id].bc_type[2][0];
            int proc_l_id   = mesh[block_id].bc_type[2][1];
            
            double**** left_states;
            left_states     = new double***[mesh[block_id].nx-1];

            for(int i=0; i<mesh[block_id].nx-1; i++){

                left_states[i]      = new double**[1];
                left_states[i][0]   = new double*[mesh[block_id].nz-1];
                    
                for(int k=0; k<mesh[block_id].nz-1; k++){
            
                    left_states[i][0][k] = new double[5];
            
                }

            }

            if(proc_l!=rank){
            
                MPI_Request* req;
                MPI_Irecv(&(left_states[0][0][0][0]), (mesh[block_id].nx-1)*(mesh[block_id].nz-1)*5, MPI_DOUBLE, proc_l, rank*1000+block_id, MPI_COMM_WORLD, req);
                MPI_Wait(req, MPI_STATUS_IGNORE);
                
            }

            for(int i=0; i<mesh[block_id].nx-1; i++){
                    
                for(int k=0; k<mesh[block_id].nz-1; k++){
            
                    if(proc_l==rank){

                        for(int sid=0; sid<5; sid++){

                            left_states[i][0][k][sid] = mesh[proc_l_id].states[i][mesh[proc_l_id].ny-2][k][sid];

                        }

                    }

                    Roe_flux( left_states[i][0][k],\
                              mesh[ block_id ].states[i][0][k],\
                              mesh[ block_id ].yn[i][0][k],\
                              mesh[ block_id ].yfluxes[i][0][k]);
            
                }

            }

            delete[] left_states;
            
        }
        
        // Sending and receiving data from top to the bottom cells
        if (mesh[block_id].bc_type[2][0]>=0){
        
            int proc_l      = mesh[block_id].bc_type[2][0];
            int proc_l_id   = mesh[block_id].bc_type[2][1];
            
            if(proc_l!=rank){

                int siz_s[4]  = {mesh[block_id].nx-1, mesh[block_id].ny-1, mesh[block_id].nz-1, 5};
                int siz_sa[4] = {mesh[block_id].nx-1, 1, mesh[block_id].nz-1, 5};
                int starts[4] = {0,0,0,0};
                MPI_Datatype Left_boundary;
                MPI_Type_create_subarray(4, siz_s, siz_sa, starts, MPI_ORDER_C, MPI_DOUBLE, &Left_boundary);
                MPI_Type_commit(&Left_boundary);
                MPI_Request* req;
                MPI_Isend(&(mesh[block_id].states[0][0][0][0]), 1, Left_boundary, proc_l, proc_l*1000+proc_l_id, MPI_COMM_WORLD, req);
            
            }

        }

        if (mesh[block_id].bc_type[3][0]>=0){

            int proc_r      = mesh[block_id].bc_type[3][0];
            int proc_r_id   = mesh[block_id].bc_type[3][1];
            
            double**** right_states;
            right_states  = new double***[mesh[block_id].ny-1];
            
            for(int i=0; i<mesh[block_id].nx-1; i++){

                right_states[i]        = new double**[1];
                right_states[i][0]     = new double*[mesh[block_id].nz-1];
                    
                for(int k=0; k<mesh[block_id].nz-1; k++){
            
                    right_states[i][0][k] = new double[5];
            
                }

            }

            if(proc_r!=rank){
            
                MPI_Request* req;
                MPI_Irecv(&(right_states[0][0][0][0]), (mesh[block_id].nx-1)*(mesh[block_id].nz-1)*5, MPI_DOUBLE, proc_r, rank*1000+block_id, MPI_COMM_WORLD, req);
                MPI_Wait(req, MPI_STATUS_IGNORE);
                
            }

            for(int i=0; i<mesh[block_id].nx-1; i++){
                    
                for(int k=0; k<mesh[block_id].nz-1; k++){
            
                    if(proc_r==rank){

                        for(int sid=0; sid<5; sid++){

                            right_states[i][0][k][sid] = mesh[proc_r_id].states[i][0][k][sid];

                        }

                    }

                    Roe_flux( mesh[ block_id ].states[i][mesh[block_id].ny-2][k],\
                              right_states[i][0][k],\
                              mesh[ block_id ].yn[i][mesh[block_id].ny-1][k],\
                              mesh[ block_id ].yfluxes[i][mesh[block_id].ny-1][k]);
            
                }

            }

            delete[] right_states;
            
        }

        // Sending and receiving data from back to the front cells
        if (mesh[block_id].bc_type[5][0]>=0){
        
            int proc_r      = mesh[block_id].bc_type[5][0];
            int proc_r_id   = mesh[block_id].bc_type[5][1];
            
            if(proc_r!=rank){

                int siz_s[4]  = {mesh[block_id].nx-1, mesh[block_id].ny-1, mesh[block_id].nz-1, 5};
                int siz_sa[4] = {mesh[block_id].nx-1, mesh[block_id].ny-1, 1, 5};
                int starts[4] = {0,0,mesh[block_id].nz-2,0};
                MPI_Datatype Right_boundary;
                MPI_Type_create_subarray(4, siz_s, siz_sa, starts, MPI_ORDER_C, MPI_DOUBLE, &Right_boundary);
                MPI_Type_commit(&Right_boundary);
                MPI_Request* req;
                MPI_Isend(&(mesh[block_id].states[0][0][0][0]), 1, Right_boundary, proc_r, proc_r*1000+proc_r_id, MPI_COMM_WORLD, req);
            
            }

        }

        if (mesh[block_id].bc_type[4][0]>=0){

            int proc_l      = mesh[block_id].bc_type[4][0];
            int proc_l_id   = mesh[block_id].bc_type[4][1];
            
            double**** left_states;
            left_states     = new double***[mesh[block_id].nx-1];
            
            for(int i=0; i<mesh[block_id].nx-1; i++){

                left_states[i] = new double**[mesh[block_id].ny-1];
                    
                for(int j=0; j<mesh[block_id].ny-1; j++){
            
                    left_states[i][j] = new double*[1];
                    left_states[i][j][0] = new double[5];
            
                }

            }

            if(proc_l!=rank){
            
                MPI_Request* req;
                MPI_Irecv(&(left_states[0][0][0][0]), (mesh[block_id].ny-1)*(mesh[block_id].nx-1)*5, MPI_DOUBLE, proc_l, rank*1000+block_id, MPI_COMM_WORLD, req);
                MPI_Wait(req, MPI_STATUS_IGNORE);
                
            }

            for(int i=0; i<mesh[block_id].nx-1; i++){
                    
                for(int j=0; j<mesh[block_id].ny-1; j++){
            
                    if(proc_l==rank){

                        for(int sid=0; sid<5; sid++){

                            left_states[i][j][0][sid] = mesh[proc_l_id].states[i][j][mesh[proc_l_id].nz-2][sid];

                        }

                    }

                    Roe_flux( left_states[i][j][0],\
                              mesh[ block_id ].states[i][j][0],\
                              mesh[ block_id ].zn[i][j][0],\
                              mesh[ block_id ].zfluxes[i][j][0]);
            
                }

            }

            delete[] left_states;
            
        }
        
        // Sending and receiving data from front to the back cells
        if (mesh[block_id].bc_type[4][0]>=0){
        
            int proc_l      = mesh[block_id].bc_type[4][0];
            int proc_l_id   = mesh[block_id].bc_type[4][1];
            
            if(proc_l!=rank){

                int siz_s[4]  = {mesh[block_id].nx-1, mesh[block_id].ny-1, mesh[block_id].nz-1, 5};
                int siz_sa[4] = {mesh[block_id].nx-1, mesh[block_id].ny-1, 1, 5};
                int starts[4] = {0,0,0,0};
                MPI_Datatype Left_boundary;
                MPI_Type_create_subarray(4, siz_s, siz_sa, starts, MPI_ORDER_C, MPI_DOUBLE, &Left_boundary);
                MPI_Type_commit(&Left_boundary);
                MPI_Request* req;
                MPI_Isend(&(mesh[block_id].states[0][0][0][0]), 1, Left_boundary, proc_l, proc_l*1000+proc_l_id, MPI_COMM_WORLD, req);
            
            }

        }

        if (mesh[block_id].bc_type[5][0]>=0){

            int proc_r      = mesh[block_id].bc_type[5][0];
            int proc_r_id   = mesh[block_id].bc_type[5][1];
            
            double**** right_states;
            right_states     = new double***[mesh[block_id].nx-1];
            
            for(int i=0; i<mesh[block_id].nx-1; i++){

                right_states[i] = new double**[mesh[block_id].ny-1];
                    
                for(int j=0; j<mesh[block_id].ny-1; j++){
            
                    right_states[i][j] = new double*[1];
                    right_states[i][j][0] = new double[5];
            
                }

            }

            if(proc_r!=rank){
            
                MPI_Request* req;
                MPI_Irecv(&(right_states[0][0][0][0]), (mesh[block_id].nx-1)*(mesh[block_id].ny-1)*5, MPI_DOUBLE, proc_r, rank*1000+block_id, MPI_COMM_WORLD, req);
                MPI_Wait(req, MPI_STATUS_IGNORE);
                
            }

            for(int i=0; i<mesh[block_id].nx-1; i++){
                    
                for(int j=0; j<mesh[block_id].ny-1; j++){
            
                    if(proc_r==rank){

                        for(int sid=0; sid<5; sid++){

                            right_states[i][j][0][sid] = mesh[proc_r_id].states[i][j][0][sid];

                        }

                    }

                    Roe_flux( mesh[ block_id ].states[i][j][mesh[block_id].nz-2],\
                              right_states[i][j][0],\
                              mesh[ block_id ].zn[i][j][mesh[block_id].nz-1],\
                              mesh[ block_id ].zfluxes[i][j][mesh[block_id].nz-1]);
            
                }

            }

            delete[] right_states;
            
        }



        // All other boundary conditions
        
        for(int bc=0; bc<6; bc++){
        
            // Inviscid wall

            if(mesh[block_id].bc_type[bc][0]==-1){

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

            if(mesh[block_id].bc_type[bc][0]==-2){

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

            if(mesh[block_id].bc_type[bc][0]==-3){

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

            if(mesh[block_id].bc_type[bc][0]==-4){

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

            if(mesh[block_id].bc_type[bc][0]==-5){

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



void gather_fluxes(int nblock, block* mesh, partition* map){

    bc_flux(nblock, mesh, map);

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



void calc_residuals(int nblock, block* mesh, partition* map){

    for(int block_id=0; block_id<nblock; block_id++){
    
        gather_fluxes(nblock,mesh,map);
        
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
