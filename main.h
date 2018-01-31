
#ifndef main_h
#define main_h

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
#include <petscvec.h>
#include <petscmat.h>
#include <petscsys.h>

/*--------------------------------------------------------*/
/* FUNCTIONS PROTOTYPES */
/*--------------------------------------------------------*/

/* FUNCTIONS TO SOLVE THE FLOW */
int compute_force_torque_fluxes(double** dudt, double** dvdt, double** domegadt, double** dTdt, double*** dCdt, double*** Ip_U, double*** Ip_V, double*** Ip_S,
                               double** U, double** V, double** T, double*** C, double** Us, double** Vs,double** Ts, double*** Cs,
                           	double* xg, double* yg, double* rp, double* Sp, double* II, 
				 double* F_drag, double* F_lift, double* Torque, double* Q_heat, double* Phi_species, int k, double t);
void get_ghosts(double** U, double** V, double** P, double** Told, double*** Cold, double CA0, double CB0);
void get_masks(double*** Ip_S, double*** Ip_U, double*** Ip_V, double** I_S, double** I_U, double** I_V,  double* xg, double* yg, double* rp, double* theta, double** coloring);
void get_Cs(double*** Cs, double*** Ip_S, double** Cp);
void get_Ts(double** Ts, double*** Ip_S, double* Tp);
void get_Us_Vs(double** Us, double** Vs, double*** Ip_U, double*** Ip_V, double** Up, double** Vp, double** Wp, double* xg, double* yg);
void get_Ustar(double** U, double** Ustar, double** V, double** A, double** Aold, double** Us, double** P, double** I_U, double ramp);
void get_Ustar_EE(double** U, double** Ustar, double** V, double** Aold, double** Us, double** P, double** I_U, double ramp);
void get_Vstar(double** V, double** Vstar, double** U, double** B, double** Bold, double** Vs, double** P, double** I_V, double ramp);
void get_Vstar_EE(double** V, double** Vstar, double** U, double** Bold, double** Vs, double** P, double** I_V, double ramp);
void old_poisson_solver(double** Ustar, double** Vstar, double **phi, double** R);
void poisson_solver(double** Ustar, double** Vstar, double **phi, int myrank, int nbproc);
void reset_phi(double** phi);
void update_flow(double** U, double** V, double** P, double** Ustar, double** Vstar, double** phi);
void update_temp_species(double** U, double** V,  double** T, double** Told, double** Ts, double** C_T, double** C_Told, double*** C, double*** Cold, double*** Cs, double*** C_C, double*** C_Cold, double** I_S, double ramp);
void update_temp_species_EE(double** U, double** V, double** T, double** Told, double** Ts, double** C_Told, double*** C, double*** Cold, double*** Cs, double*** C_Cold, double** I_S, double ramp);
void update_Xp(double* xg, double* yg, double* theta, double** Up, double** Vp, double** Wp,  int k);
void update_Cp(double** Cp, double*** Qmp, int k);
void update_Up(double** Up, double** Vp, double** Wp, double** F, double** G, double** M, int k);
void update_Tp(double* Tp, double** Qp, int k);

/* SOME HELPFUL FUNCTIONS */
void combine(char* destination, const char* path1, const char* path2);
void writeFile(FILE* file, double **data, int iStart, int iEnd, int jStart, int jEnd);
double* make1DDoubleArray(int arraySize);
double** make2DDoubleArray(int arraySizeX, int arraySizeY);
double*** make3DDoubleArray(int numberOfparticles, int arraySizeX, int arraySizeY);
int** make2DIntArray(int arraySizeX, int arraySizeY);
int*** make3DIntArray(int numberOfparticles, int arraySizeX, int arraySizeY);



/*--------------------------------------------------------*/
/* DATA */
/*--------------------------------------------------------*/
/*GLOBAL VARIABLE*/
/* DIMENSIONS */
static double Dp; 
static double d; 
static double H;
static double L;
static int ratio_L_d; 
static int ratio_d_Dp; 
static int ratio_Dp_h; 

/* PHYSICAL PARAMETERS */
static double nu;
static double rho_p, rho_f, rho_r;
static double cp, cf, cr;
static double Pr;
static double alpha_f;
static double* Df;
static double dH; // kJ/mol

/* FLOW */
static double Rp;
static double Um;
static int Np;
static int Ns;

/* TEMPERATURE PARAMETERS */
static double Tm0; 
static double Tp0;


/* GRID */
static int M;
static int N;
static int m;
static int n;
static double h;

/* TIME INTEGRATION */
static double t_move; 
static double dt;
static double dtau;

/* GAUSS_SEIDEL ALGORITHME */
double SORtol = 1e-6;
int SORitermax = 1e6;
double alpha = 1.98;


#define CLOSE_FILES \
fclose(fichier_position); \
fclose(fichier_forces); \
fclose(fichier_fluxes); \
fclose(fichier_U);\
fclose(fichier_V);\
fclose(fichier_P);\
fclose(fichier_T);\
fclose(fichier_CA);\
fclose(fichier_CB); \
fclose(fichier_data);\
fclose(fichier_mask);\
fclose(fichier_Tp);\
fclose(fichier_surface);\

#define DEFINE_FILES \
FILE* fichier_position = NULL; \
FILE* fichier_forces = NULL; \
FILE* fichier_fluxes = NULL; \
FILE* fichier_U = NULL; \
FILE* fichier_V = NULL; \
FILE* fichier_P = NULL; \
FILE* fichier_T = NULL; \
FILE* fichier_Tp = NULL; \
FILE* fichier_CA = NULL; \
FILE* fichier_CB = NULL; \
FILE* fichier_data = NULL; \
FILE* fichier_mask = NULL; \
FILE* fichier_surface = NULL; \

#define OPEN_FILES \
fichier_position = fopen("position.txt","w");\
fichier_forces = fopen("forces.txt","w");\
fichier_fluxes = fopen("fluxes.txt", "w"); \
fichier_U = fopen("U.txt","w"); \
fichier_V = fopen("V.txt","w"); \
fichier_P = fopen("P.txt", "w"); \
fichier_T = fopen("T.txt", "w"); \
fichier_Tp = fopen("Tp.txt", "w"); \
fichier_CA = fopen("CA.txt", "w"); \
fichier_CB = fopen("CB.txt", "w"); \
fichier_data = fopen("data.txt", "w"); \
fichier_mask = fopen("mask.txt", "w"); \
fichier_surface = fopen("surface.txt", "w"); \

#define OPEN_STATE \
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
state_file_particles_fluxes = fopen("state/state_particles_fluxes.txt","r+");


#define RECOVER_STATE \
ramp = 1.; \
printf("ramp = %f \n", ramp); \
fscanf(state_file, "%d", &iter_start);\
/*fscanf(state_file, "%d", &l_start);*/\
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
    fscanf(state_file_particles_fluxes, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", &Qp[k][0], &Qp[k][1], &Qp[k][2], &Qmp[k][0][0], &Qmp[k][0][1], &Qmp[k][0][2], &Qmp[k][1][0], &Qmp[k][1][1], &Qmp[k][1][2]);\
}\
fclose(state_file);\
fclose(state_file_U); fclose(state_file_Aold);\
fclose(state_file_V); fclose(state_file_Bold);\
fclose(state_file_P);\
fclose(state_file_Told); fclose(state_file_CTold);\
fclose(state_file_CAold); fclose(state_file_CCAold);\
fclose(state_file_CBold); fclose(state_file_CCBold);\
fclose(state_file_particles); fclose(state_file_particles_velocities);\
fclose(state_file_particles_forces); fclose(state_file_particles_fluxes);\

#define SAVE_STATE \
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
    fprintf(state_file_particles_fluxes, "%3.13e \t %3.13e \t %3.13e \t %3.13e \t %3.13e \t %3.13e \t %3.13e \t %3.13e \t %3.13e \n", Qp[k][0], Qp[k][1], Qp[k][2], Qmp[k][0][0], Qmp[k][0][1], Qmp[k][0][2], Qmp[k][1][0], Qmp[k][1][1], Qmp[k][1][2]);\
}\
fclose(state_file);\
fclose(state_file_U); fclose(state_file_Aold);\
fclose(state_file_V); fclose(state_file_Bold);\
fclose(state_file_P);\
fclose(state_file_Told); fclose(state_file_CTold);\
fclose(state_file_CAold); fclose(state_file_CCAold);\
fclose(state_file_CBold); fclose(state_file_CCBold);\
fclose(state_file_particles); fclose(state_file_particles_velocities);\
fclose(state_file_particles_forces); fclose(state_file_particles_fluxes);


#define WRITE_DATA \
if (fichier_data != NULL){\
fprintf(fichier_data,"d\t %f\n",d);\
fprintf(fichier_data,"H\t %f\n",H);\
fprintf(fichier_data,"ratio_L_d\t %d\n",ratio_L_d);\
fprintf(fichier_data,"nu\t %1.8f\n",nu);\
fprintf(fichier_data,"rho_p\t %f\n",rho_p);\
fprintf(fichier_data,"rho_f\t %f\n",rho_f);\
fprintf(fichier_data,"Um\t %f\n", Um); \
fprintf(fichier_data,"dt\t %f\n", dt); \
fprintf(fichier_data,"dtau\t %f\n", dtau); \
fprintf(fichier_data,"n\t %d\n",n);\
fprintf(fichier_data,"m\t %d\n",m);\
fprintf(fichier_data,"h\t %f\n",h);\
fprintf(fichier_data,"Kmax\t %d\n", Kmax);\
fprintf(fichier_data,"nKmax\t %d\n", nKmax);\
fprintf(fichier_data,"Nwrite\t %d\n", N_write);\
fprintf(fichier_data,"Twrite\t %d\n", T_write);\
fprintf(fichier_data,"Tf\t %f\n", Tf);\
fprintf(fichier_data, "Np \t %d \n", Np); \
fprintf(fichier_data, "Ns \t %d \n", Ns); \
fprintf(fichier_data, "Tm0 \t %f \n", Tm0); \
fprintf(fichier_data, "Tp0 \t %f \n", Tp0); \
fprintf(fichier_data, "Dp \t %f \n", Dp); \
fflush(fichier_data);\
}\
else{\
printf("An error occured in writeFile : invalid file pointer.\n");\
}

#endif
