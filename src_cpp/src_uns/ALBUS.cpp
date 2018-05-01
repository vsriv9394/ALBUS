#include <iostream>
#include <fstream>
#include "var_struct.cpp"
#include "read_mesh_alb.cpp"
#include "get_flux.cpp"

using namespace std;

int main(int argc, char** argv){

    double et = 0.0001;
    double dt = 10.000;

    long nodesize, facesize, cellsize, bcsize;
    char* meshname = argv[1];
    char* bcname   = argv[2];

    Node* node;
    Face* face;
    Cell* cell;
    BC*   bc;

    getmesh(nodesize, facesize, cellsize, bcsize, 
            meshname, bcname, 
            node, face, cell, bc);

    Cell* temp = new Cell[cellsize];

    ///////////////////////////////////////////////////////////////////////////

    for(double t=0.0; t<=et; t=t+dt){

        for(int i=0; i<cellsize; i++){

            temp[i].states = new double[5];

            for(int j=0; j<5; j++){

                temp[i].states[j] = cell[i].states[j];

            }

        }

        double aRK[5] = {0.000000000000000, 
                        -0.417890474499852, 
                        -1.192151694642677, 
                        -1.697784692471528, 
                        -1.514183444257156};

        double bRK[5] = {0.149659021999229,  
                         0.379210312999627, 
                         0.822955029386982, 
                         0.699450455949122, 
                         0.153057247968152};

        for(int RKiter=0; RKiter<5; RKiter++){

            for(int i=0; i<facesize; i++){

                getFlux(&face[i], temp[face[i].getLeftCellGID()].states, 
                                  temp[face[i].getRightCellGID()].states, bc);

                for(int j=0; j<5; j++){

                    temp[face[i].getLeftCellGID()].states[j] = 
                    
                        aRK[RKiter] * 
                        temp[face[i].getLeftCellGID()].states[j] -
                        face[i].flux[j] * 
                        face[i].getArea() / 
                        cell[face[i].getLeftCellGID()].getVolume() * 
                        dt;

                    cell[face[i].getLeftCellGID()].states[j] = 
                    
                        cell[face[i].getLeftCellGID()].states[j] +
                        bRK[RKiter] *
                        temp[face[i].getLeftCellGID()].states[j];



                    temp[face[i].getRightCellGID()].states[j] = 
                    
                        aRK[RKiter] * 
                        temp[face[i].getRightCellGID()].states[j] +
                        face[i].flux[j] * 
                        face[i].getArea() / 
                        cell[face[i].getRightCellGID()].getVolume() * 
                        dt;

                    cell[face[i].getRightCellGID()].states[j] = 
                    
                        cell[face[i].getRightCellGID()].states[j] +
                        bRK[RKiter] *
                        temp[face[i].getRightCellGID()].states[j];

                }

            }
        
        }

    }

    ///////////////////////////////////////////////////////////////////////////

}