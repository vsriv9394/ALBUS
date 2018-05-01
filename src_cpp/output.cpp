#include <iostream>
#include <fstream>
#include <varstruct.cpp>

using namespace std;

void writeout(char* outputfile, int nblock_local, int nblock_global, block* mesh, partition* map){

    int rank;
    int size;

    MPI_File outfile;
    MPI_Status* status = MPI_STATUS_IGNORE;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    MPI_File_open(MPI_COMM_WORLD, outputfile, MPI_MODE_CREATE || MPI_MODE_WRONLY, MPI_INFO_NULL, &outfile);

    for(int gbid=0; gbid<nblock_global; gbid++){

        for(int i=0; i<map[gbid].nbx; i++){

            for(int j=0; j<map[gbid].nby; j++){

                for(int k=0; k<map[gbid].nbz; k++){

                    ///////////////////////////////////////////////////////////////////////////////////////

                    if(map[gbid].proc[i][j][k]==rank){

                        int bid = map[gbid].proc_id[i][j][k];

                        for(int ii=0; ii<mesh[bid].nx; ii++){

                            for(int jj=0; jj<mesh[bid].ny; jj++){

                                for(int kk=0; kk<mesh[bid].nz; kk++){

                                    MPI_File_write(outfile, mesh[bid].cc[ii][jj][kk], 3, MPI_DOUBLE, status);

                                    MPI_File_write(outfile, mesh[bid].states[ii][jj][kk], 5, MPI_DOUBLE, status);
                                    MPI_File_seek(outfile, 5*sizeof(MPI_DOUBLE), MPI_SEEK_CUR);

                                }

                            }

                        }

                    }

                    else{

                        int bid = map[gbid].proc_id[i][j][k];

                        MPI_File_seek(outfile, mesh[bid].nx*mesh[bid].ny*mesh[bid].nz*8*sizeof(MPI_DOUBLE), MPI_SEEK_CUR);

                    }

                    //////////////////////////////////////////////////////////////////////////////////////

                }

            }

        }

    }

}