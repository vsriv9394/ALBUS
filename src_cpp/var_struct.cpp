using namespace std;

struct block{

    int nx,ny,nz;

    double*** x;
    double*** y;
    double*** z;

    double**** xn;
    double**** yn;
    double**** zn;

    double**** xfc;
    double**** yfc;
    double**** zfc;

    double***  xArea;
    double***  yArea;
    double***  zArea;

    double**** cc;
    double***  Volume;

    double**** xfluxes;
    double**** yfluxes;
    double**** zfluxes;
    
    double**** states;
    double**** res;
    double**** tempres;

    int bc_type[6];

    double bc_states[6][5];

};
