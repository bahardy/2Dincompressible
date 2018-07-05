
#ifndef main_h
#define main_h

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

//#define TEMP
#define MOVE
#define TWO_WAY
//#define RAMPING
#define WRITE
#define DISK
#define SLIP
//#define EXPLICIT
#define GRAVITY
#define SMOOTHING
//#define ITERATIVE
#define SEDIMENTATION

/*--------------------------------------------------------*/
/* INCLUDE HEADERS */
/*--------------------------------------------------------*/
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include <string.h>
#include <math.h>
#include <petsc.h>
#include <petscksp.h>

/*--------------------------------------------------------*/

typedef struct Data Data;
struct Data{
    /*Dimensions*/
    double Dp;
    double eps;
    double a;
    double b;
    double d;
    double H;
    double L;


    /*Physical parameters */
    double nu;
    double g;
    double rho_p, rho_f, rho_r;
    double cp, cf, cr;
    double Rep;
    double Pr;
    double Sc;
    double Le;
    double Fr;
    double Ga;
    double alpha_f;
    double* Df;
    double dH; // kJ/mol

    /*Flow */
    double u_m; // Mean streamwise velocity
    double u_max; // Max Poiseuille velocity
    int Np; // Number of particles
    int Ns; // Number of species
    double Tm0; // Initial fluid temperature
    double Tp0; // Initial particle temperature
    double CA0; // Inlet A concentration
    double CB0; // Inlet B concentration

    /*Grid */
    int M;
    int N;
    int m;
    int n;
    double h;

    /*Time integration */
    double CFL;
    double CFL_max;
    double** CFL_array;
    double r;
    double ramp;
    int Kmax;
    int nKmax;
    double t_move;
    double t_transfer;
    double dt;
    double dtau;
    double ratio_dtau_dt;
    double Tf;
    int iter;
    double Reh_max;
    double Reh_omega_max;

    /*Writings */
    int T_write;
    int N_write;

    /*Fields*/
    double** coloring;
    double*** C_n;
    double*** C_n_1;
    double* C0;
    double*** Cs;
    double** Cp;
    double* dudt;
    double* dvdt;
    double* domegadt;
    double* dTdt;
    double** dCdt;

    double** F;
    double** G;
    double** Mz;

    double* Fx;
    double* Fy;
    double* Tz;

    double** Fx_coll;
    double** Fy_coll;

    double* J;

    double*** chi_S;
    double*** chi_U;
    double*** chi_V;
    double*** Ip_S;
    double*** Ip_U;
    double*** Ip_V;
    double** I_S;
    double** I_U;
    double** I_V;


    double** P;
    double** phi;

    double** Qm;
    double*** PP;

    double* Q;
    double** QQ;
    double** Qr;
    double* rp;

    double** Reh;
    double** Reh_omega;
    double* Sp;

    double** T_n;
    double** T_n_1;
    double* Tp;
    double** Ts;


    double** u_n;
    double** u_n_1;
    double** Up;
    double** u_s;
    double** u_star;

    double** v_n;
    double** v_n_1;
    double** Vp;
    double** v_s;
    double** v_star;

    double** H_u_n_1;
    double** H_v_n_1;
    double** H_T_n_1;
    double*** H_C_n_1;

    double** omega;

    double** Omega_p;
    double** xg;
    double** yg;
    double** theta;

    double* dp;

    double **tau_xx, **tau_xy, **tau_yy;
    double Omega; 
    double R1; 
    double R2;

    double SORtol;
    int SORitermax;
    double alpha;

};


void diagnostic(Data* data);
double check_convergence(Data* data, double* Xp_old, double* Yp_old, double* theta_old,
                         double* Up_old, double* Vp_old, double* Omega_p_old);
void update_quantities(Data* data);

/* SOME HELPFUL FUNCTIONS */
void writeData(FILE* fichier_data, Data data);
double* make1DDoubleArray(int arraySize);
double** make2DDoubleArray(int arraySizeX, int arraySizeY);
double*** make3DDoubleArray(int numberOfparticles, int arraySizeX, int arraySizeY);
int** make2DIntArray(int arraySizeX, int arraySizeY);
int*** make3DIntArray(int numberOfparticles, int arraySizeX, int arraySizeY);
void free2Darray(double** array, int dim1);
void free3Darray(double*** array, int dim1, int dim2);


#endif

