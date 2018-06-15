
#ifndef main_h
#define main_h
#define TEMP
#define MOVE

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
    double* Fx;
    double* Fy;
    double** G;
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
    double** Mz;

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
    double* theta;

    double** T_n;
    double** T_n_1;
    double* Tp;
    double** Ts;

    double* Tz;

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
    double* xg;
    double* yg;
    double* dp;

    double Omega; 
    double R1; 
    double R2;

    double SORtol;
    int SORitermax;
    double alpha;

};


/* FUNCTIONS PROTOTYPES */
/*--------------------------------------------------------*/

/* FUNCTIONS TO SOLVE THE FLOW */
void compute_Qr(double** Qr, double rate, double dH, int k);
void diagnostic(Data* data);
void get_ghosts(Data* data, double T0, double* C0);
void get_masks(Data* data);
void get_Cs(Data* data);
void get_Ts(Data* data);
void get_Us_Vs(Data* data);
void get_Ustar_Vstar(Data* data, double ramp);
void get_vorticity(Data* data);
void update_flow(Data* data);
void update_scalars(Data* data);
void update_Xp(Data* data, int k);
void update_Up(Data* data, int k);
void update_Tp(Data* data,int k);
void update_Cp(Data* data, int k);
void set_up(Data* data, int argc, char *argv[], int rank);

/* SOME HELPFUL FUNCTIONS */
//void writeFile(FILE* file, double **data, int iStart, int iEnd, int jStart, int jEnd);
void writeData(FILE* fichier_data, Data data);
double* make1DDoubleArray(int arraySize);
double** make2DDoubleArray(int arraySizeX, int arraySizeY);
double*** make3DDoubleArray(int numberOfparticles, int arraySizeX, int arraySizeY);
int** make2DIntArray(int arraySizeX, int arraySizeY);
int*** make3DIntArray(int numberOfparticles, int arraySizeX, int arraySizeY);
void free2Darray(double** array, int dim1);
void free3Darray(double*** array, int dim1, int dim2);



/*#define OPEN_STATE \
FILE* state_file = NULL;\
FILE* state_file_U = NULL;\
FILE* state_file_Aold = NULL;\
FILE* state_file_V = NULL;\
FILE* state_file_Bold = NULL;\
FILE* state_file_P = NULL;\
FILE* state_file_Told = NULL;\
FILE* state_file_CTold = NULL;\
FILE* state_file_CAold = NULL;\
FILE* state_file_CCAold = NULL;\
FILE* state_file_CBold = NULL;\
FILE* state_file_CCBold = NULL;\
FILE* state_file_particles = NULL;\
FILE* state_file_particles_velocities = NULL;\
FILE* state_file_particles_forces = NULL;\
FILE* state_file_particles_fluxes = NULL;\
state_file = fopen("state/state.txt","r+");\
if(state_file == NULL){\
if(rank==0){\
printf("No state file or wrong path \n"); \
}\
}\
state_file_U = fopen("state/state_U.txt","r+");\
state_file_Aold = fopen("state/state_Aold.txt","r+");\
state_file_V = fopen("state/state_V.txt","r+");\
state_file_Bold = fopen("state/state_Bold.txt","r+");\
state_file_P = fopen("state/state_P.txt","r+");\
state_file_Told = fopen("state/state_Told.txt","r+");\
state_file_CTold = fopen("state/state_CTold.txt","r+");\
state_file_CAold = fopen("state/state_CAold.txt","r+");\
state_file_CCAold = fopen("state/state_CCAold.txt","r+");\
state_file_CBold = fopen("state/state_CBold.txt","r+");\
state_file_CCBold = fopen("state/state_CCBold.txt","r+"); \
state_file_particles = fopen("state/state_particles.txt","r+");\
state_file_particles_velocities = fopen("state/state_particles_velocities.txt","r+");\
state_file_particles_forces = fopen("state/state_particles_forces.txt","r+");\
state_file_particles_fluxes = fopen("state/state_particles_fluxes.txt","r+");*/


/*#define RECOVER_STATE \
ramp = 1.; \
printf("ramp = %f \n", ramp); \
fscanf(state_file, "%d", &iter_start);\
fscanf(state_file, "%lf", &t_start); \
\
for(int i=0; i<m; i++){\
    for(int j=0; j<n; j++){\
        fscanf(state_file_U, "%lf", &U[i][j]); fscanf(state_file_Aold, "%lf", &Aold[i][j]);\
        Ustar[i][j] = U[i][j]; \
        fscanf(state_file_V, "%lf", &V[i][j]); fscanf(state_file_Bold, "%lf", &Bold[i][j]);\
        fscanf(state_file_P, "%lf", &P[i][j]);\
        fscanf(state_file_Told, "%lf", &Told[i][j]); fscanf(state_file_CTold, "%lf", &C_Told[i][j]);\
        T[i][j] = Told[i][j]; \
        fscanf(state_file_CAold, "%lf", &Cold[0][i][j]); fscanf(state_file_CCAold, "%lf", &C_Cold[0][i][j]);\
        C[0][i][j] = Cold[0][i][j]; \
        fscanf(state_file_CBold, "%lf", &Cold[1][i][j]); fscanf(state_file_CCBold, "%lf", &C_Cold[1][i][j]);\
        C[1][i][j] = Cold[1][i][j]; \
    }\
}\
for(int k=0; k<Np; k++){\
    fscanf(state_file_particles, "%lf %lf %lf %lf %lf %lf", &xg[k], &yg[k], &rp[k], &Tp[k], &Cp[k][0], &Cp[k][1]);\
    fscanf(state_file_particles_velocities, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", &Up[k][0], &Up[k][1], &Up[k][2], &Vp[k][0], &Vp[k][1], &Vp[k][2], &Wp[k][0], &Wp[k][1], &Wp[k][2]);\
    fscanf(state_file_particles_forces, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", &F[k][0], &F[k][1], &F[k][2], &G[k][0], &G[k][1], &G[k][2], &M[k][0], &M[k][1], &M[k][2]);\
    fscanf(state_file_particles_fluxes, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", &QQ[k][0], &QQ[k][1], &QQ[k][2], &PP[k][0][0], &PP[k][0][1], &PP[k][0][2], &PP[k][1][0], &PP[k][1][1], &PP[k][1][2]);\
}\
fclose(state_file);\
fclose(state_file_U); fclose(state_file_Aold);\
fclose(state_file_V); fclose(state_file_Bold);\
fclose(state_file_P);\
fclose(state_file_Told); fclose(state_file_CTold);\
fclose(state_file_CAold); fclose(state_file_CCAold);\
fclose(state_file_CBold); fclose(state_file_CCBold);\
fclose(state_file_particles); fclose(state_file_particles_velocities);\
fclose(state_file_particles_forces); fclose(state_file_particles_fluxes);\*/

/*#define SAVE_STATE \
state_file = fopen("state/state.txt","w"); \
state_file_U = fopen("state/state_U.txt","w");\
state_file_Aold = fopen("state/state_Aold.txt","w");\
state_file_V = fopen("state/state_V.txt","w");\
state_file_Bold = fopen("state/state_Bold.txt","w");\
state_file_P = fopen("state/state_P.txt","w");\
state_file_Told = fopen("state/state_Told.txt","w");\
state_file_CTold = fopen("state/state_CTold.txt","w");\
state_file_CAold = fopen("state/state_CAold.txt","w");\
state_file_CCAold = fopen("state/state_CCAold.txt","w");\
state_file_CBold = fopen("state/state_CBold.txt","w");\
state_file_CCBold = fopen("state/state_CCBold.txt","w"); \
state_file_particles = fopen("state/state_particles.txt","w");\
state_file_particles_velocities = fopen("state/state_particles_velocities.txt","w");\
state_file_particles_forces = fopen("state/state_particles_forces.txt","w");\
state_file_particles_fluxes = fopen("state/state_particles_fluxes.txt","w"); \
\
fprintf(state_file, "%d \t %lf", iter, t); \
\
for(int i=0; i<m; i++)\
{\
    for(int j=0; j<n; j++)\
    {\
        fprintf(state_file_U, "%3.13e \t", U[i][j]); fprintf(state_file_Aold, "%3.13e \t", Aold[i][j]);\
        fprintf(state_file_V, "%3.13e \t", V[i][j]); fprintf(state_file_Bold, "%3.13e \t", Bold[i][j]);\
        fprintf(state_file_P, "%3.13e \t", P[i][j]);\
        fprintf(state_file_Told, "%3.13e \t", Told[i][j]); fprintf(state_file_CTold, "%3.13e \t", C_Told[i][j]);\
        fprintf(state_file_CAold, "%3.13e \t", Cold[0][i][j]); fprintf(state_file_CCAold, "%3.13e \t", C_Cold[0][i][j]);\
        fprintf(state_file_CBold, "%3.13e \t", Cold[1][i][j]); fprintf(state_file_CCBold, "%3.13e \t", C_Cold[1][i][j]);\
    }\
    fprintf(state_file_U, "\n"); fprintf(state_file_Aold, "\n");\
    fprintf(state_file_V, "\n"); fprintf(state_file_Bold, "\n");\
    fprintf(state_file_P, "\n");\
    fprintf(state_file_Told, "\n"); fprintf(state_file_CTold, "\n");\
    fprintf(state_file_CAold, "\n"); fprintf(state_file_CCAold, "\n");\
    fprintf(state_file_CBold, "\n"); fprintf(state_file_CCAold, "\n");\
}\
\
for(int k=0; k<Np; k++) {\
    fprintf(state_file_particles, "%3.13e \t %3.13e \t %3.13e \t %3.13e \t %3.13e \t %3.13e \n", xg[k], yg[k], rp[k], Tp[k], Cp[k][0], Cp[k][1]);\
    fprintf(state_file_particles_velocities, "%3.13e \t %3.13e \t %3.13e \t %3.13e \t %3.13e \t %3.13e \t %3.13e \t %3.13e \t %3.13e \n", Up[k][0], Up[k][1], Up[k][2], Vp[k][0], Vp[k][1], Vp[k][2], Wp[k][0], Wp[k][1], Wp[k][2]);\
    fprintf(state_file_particles_forces, "%3.13e \t %3.13e \t %3.13e \t %3.13e \t %3.13e \t %3.13e \t %3.13e \t %3.13e \t %3.13e \n", F[k][0], F[k][1], F[k][2], G[k][0], G[k][1], G[k][2], M[k][0], M[k][1], M[k][2]);\
    fprintf(state_file_particles_fluxes, "%3.13e \t %3.13e \t %3.13e \t %3.13e \t %3.13e \t %3.13e \t %3.13e \t %3.13e \t %3.13e \n", QQ[k][0], QQ[k][1], QQ[k][2], PP[k][0][0], PP[k][0][1], PP[k][0][2], PP[k][1][0], PP[k][1][1], PP[k][1][2]);\
}\
fclose(state_file);\
fclose(state_file_U); fclose(state_file_Aold);\
fclose(state_file_V); fclose(state_file_Bold);\
fclose(state_file_P);\
fclose(state_file_Told); fclose(state_file_CTold);\
fclose(state_file_CAold); fclose(state_file_CCAold);\
fclose(state_file_CBold); fclose(state_file_CCBold);\
fclose(state_file_particles); fclose(state_file_particles_velocities);\
fclose(state_file_particles_forces); fclose(state_file_particles_fluxes);*/


#endif

