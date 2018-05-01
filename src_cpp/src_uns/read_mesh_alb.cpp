#include <fstream>
#include <iostream>
#include "var_struct.cpp"

using namespace std;

void getmesh(long nodesize, long facesize, long cellsize, long bcsize, char* fname, 
             char* bcname, Node* node, Face* face, Cell* cell, BC* bc){

    ifstream infile;
    infile.open(fname, ios::binary);

    ////////////////////////////////////////////////////////////////////

    infile.read((char*)&nodesize, sizeof(nodesize));
    double x, y, z;

    node = new Node[nodesize];
    
    for(int i=0; i<nodesize; i++){

        infile.read((char*)&x, sizeof(x));
        infile.read((char*)&y, sizeof(y));
        infile.read((char*)&z, sizeof(z));
        node[i].setCoords(x,y,z);

    }

    ////////////////////////////////////////////////////////////////////

    infile.read((char*)&facesize, sizeof(facesize));
    long nd, cell_l, cell_r;
    int bctype, fctype;

    face = new Face[facesize];

    for(int i=0; i<facesize; i++){

        infile.read((char*)&fctype, sizeof(fctype));
        Node* locnode = new Node[fctype];

        for(int j=0; j<fctype; j++){

            infile.read((char*)&nd, sizeof(nd));
            locnode[j] = node[nd];

        }
        
        face[i].setFaceType(fctype);

        infile.read((char*)&cell_l, sizeof(cell_l));
        infile.read((char*)&cell_r, sizeof(cell_r));
        infile.read((char*)&bctype, sizeof(bctype));
        
        if (cell_l>=0){ face[i].setLeftCellGID(cell_l); }
        if (cell_r>=0){ face[i].setRightCellGID(cell_r); }
        face[i].setBCType(bctype);

        face[i].compute(locnode);

    }

    ////////////////////////////////////////////////////////////////////

    infile.read((char*)&cellsize, sizeof(cellsize));
    long fc;
    int cltype;
    
    cell = new Cell[cellsize];

    for(int i=0; i<cellsize; i++){

        infile.read((char*)&cltype, sizeof(cltype));
        Face* locface = new Face[cltype];

        for(int j=0; j<cltype; j++){

            infile.read((char*)&fc, sizeof(fc));
            locface[j] = face[fc];

        }

        cell[i].setCellType(cltype);

        cell[i].compute(locface);

    }

    ////////////////////////////////////////////////////////////////////

    infile.close();

    ifstream infile;
    infile.open(bcname, ios::binary);

    infile.read((char*)&bcsize, sizeof(bcsize));

    bc = new BC[bcsize];

    for(int i=0; i<bcsize; i++){

        bc[i].bcstate = new double[5];
        infile.read((char*)&bc[i].bctype, sizeof(bc[i].bctype));

        infile.read((char*)&bc[i].bcstate[0], sizeof(bc[i].bcstate[0]));
        infile.read((char*)&bc[i].bcstate[1], sizeof(bc[i].bcstate[1]));
        infile.read((char*)&bc[i].bcstate[2], sizeof(bc[i].bcstate[2]));
        infile.read((char*)&bc[i].bcstate[3], sizeof(bc[i].bcstate[3]));
        infile.read((char*)&bc[i].bcstate[4], sizeof(bc[i].bcstate[4]));

    }

    infile.close();

}