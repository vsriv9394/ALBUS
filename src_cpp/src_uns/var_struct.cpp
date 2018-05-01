// UNSTRUCTURED HEXAHEDRAL SOLVER

#include "leggauss.cpp"
#include <cmath>
#include <mpi.h>

using namespace std;

class Geometry{

    private:
        long numcells   =   0;
        long numfaces   =   0;
        long numnodes   =   0;

        int order       =   0;
        int numqp       =   1;

        int numbc       =   0;

        LegGauss lg(order, numqp);

        Node* node      =   NULL;
        Face* face      =   NULL;
        Cell* cell      =   NULL;
        BCon* bcon      =   NULL;

    public:
        void setNumCells(long dummy){numcells = dummy;}
        void setNumFaces(long dummy){numfaces = dummy;}
        void setNumNodes(long dummy){numnodes = dummy;}

        void setOrder(int* dummy){order = dummy;}
        void setNumQP(int* numqp){numqp = dummy;}

        long getNumCells(){return numcells;}
        long getNumFaces(){return numfaces;}
        long getNumNodes(){return numnodes;}

        int* getOrder(){return order;}
        int* getNumQP(){return numqp;}

};

class Node{

    private:
        double  coords[3]   =   {0.0, 0.0, 0.0};

    public:
        void setCoords(double dummy1, double dummy2, double dummy3){

            coords[0] = dummy1;
            coords[1] = dummy2;
            coords[2] = dummy3;

        }
        
        double* getCoords(){ return coords; }
        
};

class Face{

    private:
        double** left_states = NULL;
        double** left_model_states = NULL;
        double** right_states = NULL;
        double** right_model_states = NULL;

        double**  Jacobian_X = NULL;
        double**  Jacobian_Y = NULL;
        double**  Jacobian_Z = NULL;

        double** normals = NULL;

        int      bctype = 0;
        double** fluxes = NULL;
        double** model_fluxes = NULL;

    public:

        

        double** getLeftStates()       { return left_states; }
        double** getLeftModelStates()  { return left_model_states; }
        double** getRightStates()      { return right_states; }
        double** getRightModelStates() { return right_model_states; }

        double** getJacobianX()   { return JacobianX; }
        double** getJacobianY()   { return JacobianY; }
        double** getJacobianZ()   { return JacobianZ; }

        double** getNormals()     { return normals; }
        double** getFluxes()      { return fluxes; }
        double** getModelFluxes() { return model_fluxes; }

        void setLeftStates(Geometry* geometry, Equation* equation, double** leftcoeffs, double** rightcoeffs){

            for(int i=0; i<geometry->getNumQP()*geometry->getNumQP(); i++){

                for(int j=0; j<equation->getNumStates(); j++){

                    left_states[i][j] = 0;

                    for(int k=0; k<(geometry->getOrder()+1)*(geometry->getOrder()+1); k++){

                        left_states[i][j] = left_states[i][j] + leggaussLgVal

                    }

                }

            }

        }

        void    setRightCellGID(long dummy)   { right = dummy; }
        void    setBCType(int dummy)          { bctype = dummy; }
        void    setFaceType(int dummy)        { facetype = dummy; }

};

class Cell{

    private:
        double volume = 0.0;
        Node   cell_centroid;
        int    celltype = 4;

    public:
        double* states;

        void compute(Face* fc){

            double xc = 0.0;
            double yc = 0.0;
            double zc = 0.0;

            double volume_n = 0.0;

            for(int i=0; i<celltype; i++){

                xc = xc + fc[i].getFaceCentroid().getCoords()[0];    
                yc = yc + fc[i].getFaceCentroid().getCoords()[1];    
                zc = zc + fc[i].getFaceCentroid().getCoords()[2];

            }

            cell_centroid.setCoords(xc/6.0, yc/6.0, zc/6.0);

            for(int i=0; i<celltype; i++){

                volume_n = volume_n + fc[i].getArea() * abs( ( fc[i].getNormal()[0] * ( fc[i].getFaceCentroid().getCoords()[0] - xc/6.0 ) ) 
                                                           + ( fc[i].getNormal()[1] * ( fc[i].getFaceCentroid().getCoords()[1] - yc/6.0 ) ) 
                                                           + ( fc[i].getNormal()[2] * ( fc[i].getFaceCentroid().getCoords()[2] - zc/6.0 ) ) ) * ( 1/3 );

            }

            volume = volume_n;

        }

        double  getVolume()         { return volume; }
        Node*   getCellCentroid()   { return &cell_centroid; }
        int     getCellType()       { return celltype; }

        void    setCellType(int dummy)        { celltype = dummy; }
        
};

class BC{

    public:
        double* bcstate;
        int bctype;

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class Equation{

    private:
        int     num_states              = 1;
        int     num_model_states        = 0;
        double  dt;
        double  et;

    public:
        void (*Time_Int) (double, Cell*, Cell*);

        int*    getOrder()          {return order;}
        int*    getNumQuadPoints()  {return num_quad_points;}
        double  getDt()             {return dt;}
        double  getEt()             {return et;}
        int     getNumStates()      {return num_states;}
        int     getNumModelStates() {return num_model_states;}

        void setOrder           (int* dummy)    {order = dummy;}
        void setNumQuadPoints   (int* dummy)    {num_quad_points = dummy;}
        void setDt              (double dummy)  {dt = dummy;}
        void setEt              (double dummy)  {et = dummy;}

};

class Euler : public Equation{

    states = 5;
    
    map<string,void(*)(int,Face*,Cell*)> Eulerflux;

    EulerFlux["Roe"]        = Roe_flux;
    EulerFlux["HLLC"]       = HLLC;

    public:
        void updateEuler(string scheme, int numFace, Face* face, Cell* cell){
            
            EulerFlux[scheme](numFace,face,cell);
            
        }

};