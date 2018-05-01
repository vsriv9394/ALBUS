#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "var_struct.cpp"

using namespace std;

void read_mesh(char* meshfile, int& nblock_local, int& nblock_global, block* mesh, partition* map){

    // Initialize number of blocks on this processor to zero
    
    nblock_local = 0;

    int rank;   // Rank of this processor
    int size;   // Number of processors in the MPI proc range

    MPI_File mfile;                         // File variable for the mesh file
    MPI_Status* status = MPI_STATUS_IGNORE; // Dummy MPI status
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);   // Assign the rank variable
    MPI_Comm_size(MPI_COMM_WORLD, &size);   // Assign the size variable

    // Opening the mesh file
    
    MPI_File_open(MPI_COMM_WORLD, meshfile, MPI_MODE_RDONLY, MPI_INFO_NULL, &mfile);
    
    // int nblock_global;                                      // Number of blocks in the mesh file
    MPI_File_read(mfile,&nblock_global,1,MPI_INT,status);   // Read "nblock_global" from the mesh file
    
    map = new partition[nblock_global];                     // Initialize partition array "map"
    
    int maxprocs;                                           // Total processors required for the mesh after partitioning
    MPI_File_read(mfile,&maxprocs,1,MPI_INT,status);        // Read "maxprocs" from the mesh file

    // Initialize the number of blocks currently on each of the "maxprocs" number of processors to zero
    
    int* proc_local_ind = new int[maxprocs];

    for(int procid=0; procid<maxprocs; procid++){

        proc_local_ind[procid] = 0;

    }

    // Display a notification if the available processing resources are more than required

    if(rank==0 && size>maxprocs){
        
        cout<<"Number of processors available exceeds the required number!!"<<endl;
    
    }

    // Display a notification if the available processing resources are less than required

    if(rank==0 && size<maxprocs){

        cout<<"Number of processors available short of the required number!!"<<endl;
        return;

    }

    // If the available resources are sufficient then execute the else condition
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    else{

        vector<int> global_block_id;        // Vector containing indices of global mesh blocks corresponding to the local blocks on this processor
        vector<int> ibx;                    // Vector containing x-indices of local processor blocks which are sub-blocks relative to the global mesh blocks with indices "global_block_id" 
        vector<int> iby;                    // Vector containing y-indices of local processor blocks which are sub-blocks relative to the global mesh blocks with indices "global_block_id" 
        vector<int> ibz;                    // Vector containing z-indices of local processor blocks which are sub-blocks relative to the global mesh blocks with indices "global_block_id" 

        /////////////////////////////////////////
        // Setup processor-block mapping
        /////////////////////////////////////////

        // Loop over all the global blocks in the mesh file

        for(int m=0; m<nblock_global; m++){

            int proci, procd;

            MPI_File_read(mfile,&proci,1,MPI_INT,status);           // Read rank of processor for sub-block [0,0,0]
            MPI_File_read(mfile,&procd,1,MPI_INT,status);           // Read the jump in processor rank for subsequent sub-blocks
            MPI_File_read(mfile,&map[m].nbx,1,MPI_INT,status);      // Number of sub-blocks in the x-direction
            MPI_File_read(mfile,&map[m].nby,1,MPI_INT,status);      // Number of sub-blocks in the y-direction
            MPI_File_read(mfile,&map[m].nbz,1,MPI_INT,status);      // Number of sub-blocks in the z-direction

            // Setup 3D-array for processor rank ("proc") and local index of the current sub-block on that processor ("proc_id") 
            
            map[m].proc     = new int**[map[m].nbx];
            map[m].proc_id  = new int**[map[m].nbx];
            
            for(int i=0; i<map[m].nbx; i++){

                map[m].proc[i]     = new int*[map[m].nby];
                map[m].proc_id[i]  = new int*[map[m].nby];

                for(int j=0; j<map[m].nby; j++){

                    map[m].proc[i][j]     = new int[map[m].nbz];
                    map[m].proc_id[i][j]  = new int[map[m].nbz];

                }

            }

            // Assign values to the processor-block mapping

            for(int k=0; k<map[m].nbz; k++){

                for(int j=0; j<map[m].nby; j++){

                    for(int i=0; i<map[m].nbx; i++){

                        map[m].proc[i][j][k]    = proci;                    // Assign the processor rank for the current sub-block
                        map[m].proc_id[i][j][k] = proc_local_ind[proci];    // Assign the number of local blocks already existing in the processor as index
                        proc_local_ind[proci]++;                            // Increment the number of local blocks in the processor

                        if(proci==rank){

                            nblock_local++;                                 // Increment the number of local blocks in this processor
                            global_block_id.push_back(m);                   // Push a new value for "global_block_id" as "m"
                            ibx.push_back(i);                               // Push a new value for "ibx" as iterate "i"
                            iby.push_back(j);                               // Push a new value for "iby" as iterate "j"
                            ibz.push_back(k);                               // Push a new value for "ibz" as iterate "k"

                        }

                        proci = proci+procd;                                // Make the jump in processor rank as required for the next sub-block

                    }

                }

            }

        }

        delete[] proc_local_ind;                                            // No use of the number of local blocks on each processor now

        ///////////////////////////////////////////////////////////

        // Setup local block array

        block* mesh = new block[nblock_local];

        ///////////////////////////////////////////////////////////

        // Set the boundary conditions for local blocks

        int count = 0;                                                      // ***IMPORTANT***: count defines the "index" of local blocks on this processor

        // Looping over global mesh blocks

        for(int m=0; m<nblock_global; m++){

            // Checking if a sub-block of this global block is among the local blocks on this processor
            
            if(global_block_id[count]==m){

                MPI_File_read(mfile,mesh[count].bc_type[0],  6, MPI_INT   ,status);     // Read the boundary conditions for the global block
                MPI_File_read(mfile,mesh[count].bc_states,30,MPI_DOUBLE,status);        // Read the state variables required for all 6 boundaries
                MPI_File_seek(mfile,30*sizeof(MPI_DOUBLE),MPI_SEEK_CUR);                // Advance the MPI_File_read cursor to the next block's boundary conditions
                
                // Setup the boundary conditions on the local block
                ///////////////////////////   X   /////////////////////////////////////

                if(ibx[count]==0){                                                                                  // Check if it is the left-most block

                    if(mesh[count].bc_type[0][0]>=0){

                        int iblk = mesh[count].bc_type[0][0];                                                       // Which global block is adjacent to the current block
                        mesh[count].bc_type[0][0] = map[iblk].proc[map[iblk].nx-1][iby[count]][ibz[count]];         // Which processor to look for
                        mesh[count].bc_type[1][0] = map[iblk].proc_id[map[iblk].nx-1][iby[count]][ibz[count]];      // Which local block on that proc to look for

                    }

                }
                else{

                    mesh[count].bc_type[0][0] = map[m].proc[ibx[count]-1][iby[count]][ibz[count]];
                    mesh[count].bc_type[1][0] = map[m].proc_id[ibx[count]-1][iby[count]][ibz[count]];


                }

                if(ibx[count]==map[m].nbx-1){

                    if(mesh[count].bc_type[0][1]>=0){

                        int iblk = mesh[count].bc_type[0][1];
                        mesh[count].bc_type[0][1] = map[iblk].proc[0][iby[count]][ibz[count]];
                        mesh[count].bc_type[1][1] = map[iblk].proc_id[0][iby[count]][ibz[count]];

                    }

                }
                else{

                    mesh[count].bc_type[0][1] = map[m].proc[ibx[count]+1][iby[count]][ibz[count]];
                    mesh[count].bc_type[1][1] = map[m].proc_id[ibx[count]+1][iby[count]][ibz[count]];


                }

                ///////////////////////////   Y   /////////////////////////////////////

                if(iby[count]==0){

                    if(mesh[count].bc_type[0][2]>=0){

                        int iblk = mesh[count].bc_type[0][2];
                        mesh[count].bc_type[0][2] = map[iblk].proc[ibx[count]][map[iblk].ny-1][ibz[count]];
                        mesh[count].bc_type[1][2] = map[iblk].proc_id[ibx[count]][map[iblk].ny-1][ibz[count]];

                    }

                }
                else{

                    mesh[count].bc_type[0][2] = map[m].proc[ibx[count]][iby[count]-1][ibz[count]];
                    mesh[count].bc_type[1][2] = map[m].proc_id[ibx[count]][iby[count]-1][ibz[count]];


                }

                if(iby[count]==map[m].nby-1){

                    if(mesh[count].bc_type[0][3]>=0){

                        int iblk = mesh[count].bc_type[0][3];
                        mesh[count].bc_type[0][3] = map[iblk].proc[ibx[count]][0][ibz[count]];
                        mesh[count].bc_type[1][3] = map[iblk].proc_id[ibx[count]][0][ibz[count]];

                    }

                }
                else{

                    mesh[count].bc_type[0][3] = map[m].proc[ibx[count]][iby[count]+1][ibz[count]];
                    mesh[count].bc_type[1][3] = map[m].proc_id[ibx[count]][iby[count]+1][ibz[count]];

                }

                ///////////////////////////   Z   /////////////////////////////////////

                if(ibz[count]==0){

                    if(mesh[count].bc_type[0][4]>=0){

                        int iblk = mesh[count].bc_type[0][4];
                        mesh[count].bc_type[0][4] = map[iblk].proc[ibx[count]][iby[count]][map[iblk].nz-1];
                        mesh[count].bc_type[1][4] = map[iblk].proc_id[ibx[count]][iby[count]][map[iblk].nz-1];

                    }

                }
                else{

                    mesh[count].bc_type[0][4] = map[m].proc[ibx[count]][iby[count]][ibz[count]-1];
                    mesh[count].bc_type[1][4] = map[m].proc_id[ibx[count]][iby[count]][ibz[count]-1];


                }

                if(ibz[count]==map[m].nbz-1){

                    if(mesh[count].bc_type[0][5]>=0){

                        int iblk = mesh[count].bc_type[0][5];
                        mesh[count].bc_type[0][5] = map[iblk].proc[ibx[count]][iby[count]][0];
                        mesh[count].bc_type[1][5] = map[iblk].proc_id[ibx[count]][iby[count]][0];

                    }

                }
                else{

                    mesh[count].bc_type[0][5] = map[m].proc[ibx[count]][iby[count]][ibz[count]+1];
                    mesh[count].bc_type[1][5] = map[m].proc_id[ibx[count]][iby[count]][ibz[count]+1];


                }

                count++;            // Advance the local block index

            }

            else{                   // Skip this global block's boundary condition

                MPI_File_seek(mfile, 6*sizeof(MPI_INT),   MPI_SEEK_CUR);
                MPI_File_seek(mfile,30*sizeof(MPI_DOUBLE),MPI_SEEK_CUR);

            }

        }

        // Set the dimensions of block

        count = 0;

        for(int m=0; m<nblock_global; m++){

            MPI_File_read(mfile,&map[m].nx,1,MPI_INT,status);
            MPI_File_read(mfile,&map[m].ny,1,MPI_INT,status);
            MPI_File_read(mfile,&map[m].nz,1,MPI_INT,status);
                
            if(global_block_id[count]==m){

                mesh[count].nx = ((map[m].nx-1)/map[m].nbx)+1;
                mesh[count].ny = ((map[m].ny-1)/map[m].nby)+1;
                mesh[count].nz = ((map[m].nz-1)/map[m].nbz)+1;
                count++;

            }

        }

        MPI_File_seek(mfile,3*sizeof(MPI_INT),MPI_SEEK_CUR);

        // Initialize and set the coordinates of block

        count = 0;

        for(int m=0; m<nblock_global; m++){

            if(global_block_id[count]==m){

                mesh[count].x = new double**[mesh[count].nx];
                mesh[count].y = new double**[mesh[count].nx];
                mesh[count].z = new double**[mesh[count].nx];

                for(int i=0; i<mesh[count].nx; i++){

                    mesh[count].x[i] = new double*[mesh[count].ny];
                    mesh[count].y[i] = new double*[mesh[count].ny];
                    mesh[count].z[i] = new double*[mesh[count].ny];

                    for(int j=0; j<mesh[count].ny; j++){

                        mesh[count].x[i][j] = new double[mesh[count].nz];
                        mesh[count].y[i][j] = new double[mesh[count].nz];
                        mesh[count].z[i][j] = new double[mesh[count].nz];

                    }

                }

                int x1,xe,y1,ye,z1,ze;
                
                x1 = ibx[count]*(map[m].nx-1)/map[m].nbx;
                xe = (ibx[count]+1)*(map[m].nx-1)/map[m].nbx;

                y1 = iby[count]*(map[m].ny-1)/map[m].nby;
                ye = (iby[count]+1)*(map[m].ny-1)/map[m].nby;

                z1 = ibz[count]*(map[m].nz-1)/map[m].nbz;
                ze = (ibz[count]+1)*(map[m].nz-1)/map[m].nbz;

                MPI_File_seek(mfile, (z1*map[m].nx*map[m].ny + y1*map[m].nx + x1)*3*sizeof(MPI_DOUBLE), MPI_SEEK_CUR);
                
                for(int k=0; k<mesh[count].nz; k++){

                    for(int j=0; j<mesh[count].ny; j++){

                        for(int i=0; i<mesh[count].nx; i++){

                            MPI_File_read(mfile,&mesh[count].x[i][j][k],1,MPI_DOUBLE,status);
                            MPI_File_read(mfile,&mesh[count].y[i][j][k],1,MPI_DOUBLE,status);
                            MPI_File_read(mfile,&mesh[count].z[i][j][k],1,MPI_DOUBLE,status);
                            MPI_File_seek(mfile,sizeof(MPI_DOUBLE),MPI_SEEK_CUR);

                        }

                        MPI_File_seek(mfile,(map[m].nx-(xe-x1+1))*3*sizeof(MPI_DOUBLE),MPI_SEEK_CUR);

                    }

                    MPI_File_seek(mfile,(map[m].nx-(xe-x1+1) + map[m].nx*(map[m].ny-(ye-y1+1)))*3*sizeof(MPI_DOUBLE),MPI_SEEK_CUR);

                }
                
                count++;

                MPI_File_seek(mfile,(map[m].nx-xe + map[m].nx*(map[m].ny-ye-1) + map[m].ny*map[m].nx*(map[m].nz-ze-1))*3*sizeof(MPI_DOUBLE),MPI_SEEK_CUR);

            }

            else{

                MPI_File_seek(mfile,map[m].nx*map[m].ny*map[m].nz*3*sizeof(MPI_DOUBLE),MPI_SEEK_CUR);
                
                count++;

            }

        }

        // Evaluation of geometrical quantities
        
        for(int block_id=0; block_id<nblock_local; block_id++){

            cout<<"Initializing block geometrical quantities"<<block_id+1<<endl;

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

            double s1x, s2x, s3x, s4x,
                   s1y, s2y, s3y, s4y,
                   s1z, s2z, s3z, s4z;

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

    }

    MPI_File_close(&mfile);
    MPI_Barrier(MPI_COMM_WORLD);
    return;

}
