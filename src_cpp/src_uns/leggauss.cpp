using namespace std;

class LegGauss{

    private:

        double Legendre_coeff[5][5] =   
        {   {   1.0,    0.0,    0.0,    0.0,    0.0},
            {   0.0,    1.0,    0.0,    0.0,    0.0},
            {  -0.5,    0.0,    1.5,    0.0,    0.0},
            {   0.0,   -1.5,    0.0,    2.5,    0.0},
            { 0.375,    0.0,  -3.75,    0.0,  4.375} };

        double Legendre_deriv[5][5] =   
        {   {   0.0,    0.0,    0.0,    0.0,    0.0},
            {   1.0,    0.0,    0.0,    0.0,    0.0},
            {   0.0,    3.0,    0.0,    0.0,    0.0},
            {  -3.0,    0.0,    5.0,    0.0,    0.0},
            {   0.0,   -7.5,    0.0,   8.75,    0.0} };

        double LG_Quad_points[5][3] =   
        {   {           0.0,            0.0,            0.0},
            {  0.5773502692,            0.0,            0.0},
            {           0.0,  0.77459666924,            0.0},
            {  0.3399810436,  0.86113631159,            0.0},
            {           0.0,  0.53846931011,  0.90617984594} };

        double LG_Quad_weight[5][3] =   
        {   {           2.0,            0.0,            0.0},
            {           1.0,            0.0,            0.0},
            {  0.8888888889,  0.55555555556,            0.0},
            {  0.6521451549,  0.34785484514,            0.0},
            {  0.5688888889,  0.47862867050,  0.23692688506} };

    public:
        double *qp, *qc, **lv, **ld, *ev_tail, *ev_head, *ed_head, *ed_tail;

        LegGauss(int order, int numQP){

            qp      =   new double[numQP];            qc      =   new double[numQP];

            lv      =   new double*[order];           ld      =   new double*[order];
            
            ev_tail =   new double[order];            ev_head =   new double[order];
            ed_tail =   new double[order];            ed_head =   new double[order];

            for(int i=0; i<numQP; i++){

                qp[i]   =   LG_Quad_points[(numQP-i)/2] * pow(-1.0,(numQP-i)%2);
                qc[i]   =   LG_Quad_weight[(numQP-i)/2];

            }
            
            for(int i=0; i<=order; i++){

                for(int j=0; j<numQP; j++){

                    lv[i][j] = 0.0;
                    ld[i][j] = 0.0;

                    if(j==0){

                        ev_head[i] = 0.0;
                        ev_tail[i] = 0.0;
                        ed_head[i] = 0.0;
                        ed_tail[i] = 0.0;

                    }

                    for(int k=0; k<=i; k++){

                        lv[i][j] = lv[i][j] + Legendre_coeff[i][k]*pow(qp[j],k);
                        ld[i][j] = ld[i][j] + Legendre_deriv[i][k]*pow(qp[j],k);

                        if(j==0){

                            ev_head[i] = ev_head[i] + Legendre_coeff[i][k];
                            ev_tail[i] = ev_tail[i] + Legendre_coeff[i][k]*pow(-1.0,k);
                            ed_head[i] = ed_head[i] + Legendre_deriv[i][k];
                            ed_tail[i] = ed_tail[i] + Legendre_deriv[i][k]*pow(-1.0,k);

                        }

                    }

                }

            }

        }            



};