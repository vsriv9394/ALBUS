#include <iostream>
#include <fstream>
#include "calc_res.cpp"

void LSRK4(int nblock, block* mesh, double dt){

    double aRK[5] = {0.000000000000000, -0.417890474499852, -1.192151694642677, -1.697784692471528, -1.514183444257156};
    double bRK[5] = {0.149659021999229,  0.379210312999627,  0.822955029386982,  0.699450455949122,  0.153057247968152};

    for(int RKiter=0; RKiter<5; RKiter++){
    
        calc_residuals(nblock, mesh);

        for(int block_id=0; block_id<nblock; block_id++){
        
            for(int i=0; i<mesh[block_id].nx-1; i++){
            
                for(int j=0; j<mesh[block_id].ny-1; j++){
                
                    for(int k=0; k<mesh[block_id].nz-1; k++){

                        for(int l=0; l<5; l++){
                    
                            mesh[block_id].tempres[i][j][k][l] = aRK[RKiter]*\
                                                                 mesh[block_id].tempres[i][j][k][l]+\
                                                                 mesh[block_id].res[i][j][k][l]*dt;

                            mesh[block_id].states[i][j][k][l] = mesh[block_id].states[i][j][k][l]+\
                                                                bRK[RKiter]*\
                                                                mesh[block_id].tempres[i][j][k][l];

                        }

                    }
                
                }
            
            }
        
        }
    
    }

}
