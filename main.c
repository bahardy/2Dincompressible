#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <petsc.h>
#include <sys/stat.h>

#include "main.h"
#include "poisson.h"
#include "write.h"
#include "fields_memory.h"
#include "forces.h"
#include "collision.h"
#include "particle_motion.h"
#include "set_up.h"
#include "flow_solver.h"


int main(int argc, char *argv[]){

    PetscInitialize(&argc, &argv, 0, 0);
    int rank, nbproc;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank );
    MPI_Comm_size(PETSC_COMM_WORLD, &nbproc);
    printf("Hello from rank %d \n", rank);

    struct stat st = {0};
    if (stat("results", &st) == -1) {
        mkdir("results", 0700);
    }


    /**------------------------------- DATA BASE CREATION ------------------------------- **/
    Data data;
    set_up(&data, argc, argv, rank);
    FILE* fichier_data = fopen("results/data.txt", "w");
    writeData(fichier_data, data);
    fclose(fichier_data);

    /** -------- Some files creation and data writing-------- **/

    FILE** fichier_particles = malloc(sizeof(FILE*)*data.Np);
    FILE** fichier_forces = malloc(sizeof(FILE*)*data.Np);
    FILE** fichier_fluxes = malloc(sizeof(FILE*)*data.Np);
    FILE* fichier_stat = fopen("results/stats.txt", "w");
    FILE* fichier_forces_NOCA = fopen("results/forces_NOCA.txt", "w+");

    for(int k = 0; k<data.Np; k++)
    {
        char K[10];
        sprintf(K, "%d", k);

        char fileParticle[30];
        strcpy(fileParticle, "results/particle");
        strcat(fileParticle, "-");
        strcat(fileParticle, K);
        strcat(fileParticle, ".txt");
        fichier_particles[k] = fopen(fileParticle, "w+");

        char fileForces[30];
        strcpy(fileForces, "results/forces");
        strcat(fileForces, "-");
        strcat(fileForces, K);
        strcat(fileForces, ".txt");
        fichier_forces[k] = fopen(fileForces, "w+");

        char fileFluxes[30];
        strcpy(fileFluxes, "results/fluxes");
        strcat(fileFluxes, "-");
        strcat(fileFluxes, K);
        strcat(fileFluxes, ".txt");
        fichier_fluxes[k] = fopen(fileFluxes, "w+");
    }

    /** ------------------------------- FIELDS INITIALIZATION ------------------------------- **/
    allocate_fields(&data);
    initialize_fields(&data);
    get_masks(&data);
    get_Us_Vs(&data);

    /** ----- BOUNDARY CONDITION -----------**/
    get_ghosts(&data, data.Tm0, data.C0);


    /** ------------------------------- RAMPING ------------------------------- **/

#ifdef RAMPING
    for(int K=0; K<=data.nKmax*data.Kmax; K++) {
        /* role de nKmax : atteindre une solution stable avec les particules fixées avant de les "lacher" */
        data.ramp = fmin(1., (double) K / data.Kmax);

        PetscPrintf(PETSC_COMM_WORLD, "\n \n BEGIN ramp = %f \n", data.ramp);
        get_masks(&data);
        get_Ts(&data);
        get_Cs(&data);
        for (int k = 0; k < Np; k++) {
            integrate_penalization(&data, &surf, k);
            /* dudt, dvdt, etc = 0 because the particle is fixed */
            compute_forces_fluxes(&data, k);
        }

        get_Ustar_Vstar(&data, data.ramp);
        poisson_solver(&data, rank, nbproc);
        update_flow(&data);
        get_ghosts(&data, data.Tm0, data.C0);
        diagnostic(&data);

    }
#endif

#ifdef WRITE
    /*INITIAL SOLUTION (t=0) AFTER RAMPING */
    if(rank==0){
        writeFields(&data, 0);
    }
#endif

    /** -------------------------------TIME STEPPING ------------------------------- **/
    double t = 0;
    double surf = 0.;
    int k = 0;

    data.ramp = 1;
    data.iter = 1;

    // we feed reactants at the inlet
    data.C0[0] = data.CA0;
    //data.C0[1] = data.CB0;

    while(t < data.Tf){

        PetscPrintf(PETSC_COMM_WORLD, "\n \n BEGIN iter %d : t = %f \n", data.iter, t);
        for (k = 0; k<data.Np; k++) {
            PetscPrintf(PETSC_COMM_WORLD, "Position of the center of mass of particle %d: (x,y) = (%f,%f) \n", k + 1,
                        data.xg[k], data.yg[k]);
            PetscPrintf(PETSC_COMM_WORLD, "Angle: theta  = %f \n", data.theta[k]);
        }

        /** --- SOLVE SOLID PHASE --- */
        int flag_out = 0;

        /** Check for collisions **/
        collision(&data);

        for (k = 0; k<data.Np; k++){
            /* Integrate penalization term */
            flag_out += integrate_penalization(&data, &surf, k);
#ifdef  MOVE
            /* Velocity - Forces */
            if(t > data.t_move){
                update_Xp(&data, k);
                update_Up(&data, k);
            }
#endif
#ifdef  TEMP
            /*Temperature - Species - Fluxes */
            if(t > data.t_transfer)
            {
                update_Tp(&data, k);
                update_Cp(&data, k);
            }
#endif
            compute_forces_fluxes(&data, k);
        }
        if (flag_out > 0){
            break;
        }

        /** --- SOLVE FLUID PHASE --- */
#ifdef  MOVE
        /*Compute the mask functions */
        get_masks(&data);
        get_Us_Vs(&data);
#endif
#ifdef TEMP
        get_Ts(&data);
        //get_Cs(&data);
#endif
        get_Ustar_Vstar(&data, data.ramp);

        clock_t t_init = clock();
        poisson_solver(&data, rank, nbproc);
        clock_t t_final = clock();
        double t_Poisson = ((double) (t_final - t_init))/CLOCKS_PER_SEC;
        PetscPrintf(PETSC_COMM_WORLD, "Poisson solver took %f seconds \n", t_Poisson);
        poisson_residual(&data);

#ifdef TEMP
        update_scalars(&data);
#endif
        update_flow(&data);
        get_ghosts(&data, data.Tm0, data.C0);
        get_vorticity(&data);
        get_tau(&data);

        compute_forces_NOCA(&data, fichier_forces_NOCA, 1, data.m-2, 1, data.n-2);
        diagnostic(&data);
        update_quantities(&data); 


#ifdef WRITE
        if(rank == 0){
            writeStatistics(&data, fichier_stat);
            for (k = 0; k< data.Np; k++) {
                writeForces(&data, fichier_forces, k);
                writeParticle(&data, fichier_particles, k);
                writeFluxes(&data, fichier_fluxes, k);
            }
            if(data.iter % data.T_write == 0){
                writeFields(&data, data.iter);
            }
        }
#endif
        /* Increment time */
        t += data.dt;
        data.iter ++;

    }
    for (k = 0; k< data.Np; k++)
    {
        fclose(fichier_forces[k]);
        fclose(fichier_fluxes[k]);
        fclose(fichier_particles[k]);
    }
    free(fichier_forces);
    free(fichier_fluxes);
    free(fichier_particles);

    fclose(fichier_stat);
    fclose(fichier_forces_NOCA);

    free_fields(&data);

    PetscFinalize();
    return 0;
}


void diagnostic(Data* data)
{
    double** u_n = data->u_n;
    double** v_n = data->v_n;
    double** omega = data->omega;
    double** Reh = data->Reh;
    double** Reh_omega = data->Reh_omega;
    double h = data->h;
    double dt = data->dt;
    double** u = data->u_n;
    double** v = data->v_n;
    double** CFL_array = data->CFL_array;
    double nu = data->nu;
    double m = data->m;
    double n = data->n;

    data->Reh_max = 0.;
    data->Reh_omega_max = 0.;
    data->CFL_max = 0.;

    int i, j;
    for (i=1; i<m-1; i++) {
        for (j = 1; j<n-1; j++) {
            Reh[i][j] = (fabs(u_n[i][j]) + fabs(v_n[i][j]))*h/nu;
            Reh_omega[i][j] = fabs(omega[i][j])*h*h/n;
            CFL_array[i][j] = (fabs(u[i][j]) + fabs(v[i][j]))*dt/h;
            if(Reh[i][j] > data->Reh_max){
                data->Reh_max = Reh[i][j];
            }
            if(Reh_omega[i][j] > data->Reh_omega_max){
                data->Reh_omega_max = Reh_omega[i][j];
            }
            if(CFL_array[i][j] > data->CFL_max){
                data->CFL_max = CFL_array[i][j];
            }
        }
    }
}

void update_quantities(Data* data)
{
    for (int k=0; k<data->Np; k++) {

        data->Fx_coll[k][0] =  data->Fx_coll[k][1];
        data->Fx_coll[k][1] =  data->Fx_coll[k][2];
        data->Fx_coll[k][2] = 0;

        data->Fy_coll[k][0] =  data->Fy_coll[k][1];
        data->Fy_coll[k][1] =  data->Fy_coll[k][2];
        data->Fy_coll[k][2] = 0;

        /* Force along x-direction */
        data->F[k][0] = data->F[k][1]; /* n-2*/
        data->F[k][1] = data->F[k][2]; /* n-1*/

        /* Force along y-direction */
        data->G[k][0] = data->G[k][1];
        data->G[k][1] = data->G[k][2];

        /* Moment along z-direction */
        data->Mz[k][0] = data->Mz[k][1];
        data->Mz[k][1] = data->Mz[k][2];

        data->xg[k][0] = data->xg[k][1];
        data->yg[k][0] = data->yg[k][1];
        data->theta[k][0] = data->theta[k][1];

#ifdef TWO_WAY
        data->Up[k][0] = data->Up[k][1];
        data->Up[k][1] = data->Up[k][2];

        data->Vp[k][0] = data->Vp[k][1];
        data->Up[k][1] = data->Up[k][2];

        data->Omega_p[k][0] = data->Omega_p[k][1];
        data->Omega_p[k][1] = data->Omega_p[k][2];

#endif
    }
}

/*
void update_Xp(Data* data, int k)
{
    double* xg = data->xg;
    double* yg = data->yg;
    double* theta = data->theta;
    double** Up = data->Up;
    double** Vp = data->Vp;
    double** Omega_p = data->Omega_p;

    double dt = data->dt;

    xg[k] += dt*(23.*Up[k][2]-16.*Up[k][1]+5.*Up[k][0])/12.;
    yg[k] += dt*(23.*Vp[k][2]-16.*Vp[k][1]+5.*Vp[k][0])/12.;
    theta[k] += dt*(23.*Omega_p[k][2]-16.*Omega_p[k][1]+5.*Omega_p[k][0])/12.;
    PetscPrintf(PETSC_COMM_WORLD,"Position of the center of mass of particle %d: (x,y) = (%f,%f) \n", k+1, xg[k], yg[k]);
    PetscPrintf(PETSC_COMM_WORLD,"Angle: theta  = %f \n", theta[k]);

}

void update_Up(Data* data, int k)
{
    double* dudt = data->dudt;
    double* dvdt = data->dvdt;
    double* domegadt = data->domegadt;
    double rho_r = data->rho_r;
    double rho_f = data->rho_f;
    double rho_p = data->rho_p;
    double* Sp = data->Sp;
    double* J = data->J; //m^4
    double** F = data->F;
    double** G = data->G;
    double** M = data->Mz;
    double g = data->g;
    double** Up = data->Up;
    double** Vp = data->Vp;
    double** Omega_p = data->Omega_p;
    double** Fx_coll = data->Fx_coll;
    double** Fy_coll = data->Fy_coll;

    double dt = data->dt;


    dudt[k] = (23.*F[k][2]-16.*F[k][1]+5.*F[k][0])/(12.*Sp[k]*(rho_r - 1.)) + (23.*Fx_coll[k][2]-16.*Fx_coll[k][1]+5.*Fx_coll[k][0])/(12.*Sp[k]*(rho_p - rho_f)) - g;
    Up[k][3] = Up[k][2] + dt*dudt[k];
    dvdt[k] = (23.*G[k][2]-16.*G[k][1]+5.*G[k][0])/(12.*Sp[k]*(rho_r - 1.)) + (23.*Fy_coll[k][2]-16.*Fy_coll[k][1]+5.*Fy_coll[k][0])/(12.*Sp[k]*(rho_p - rho_f)) ;
    Vp[k][3] = Vp[k][2] + dt*dvdt[k];
    domegadt[k] = (23.*M[k][2]-16.*M[k][1]+5.*M[k][0])/(12.*J[k]*(rho_r - 1.));
    Omega_p[k][3] = Omega_p[k][2] + dt*domegadt[k];

    PetscPrintf(PETSC_COMM_WORLD, "Up = %f \n", Up[k][3]);
    Up[k][0] = Up[k][1]; Up[k][1] = Up[k][2]; Up[k][2] = Up[k][3];
    Vp[k][0] = Vp[k][1]; Vp[k][1] = Vp[k][2]; Vp[k][2] = Vp[k][3];
    Omega_p[k][0] = Omega_p[k][1]; Omega_p[k][1] = Omega_p[k][2]; Omega_p[k][2] = Omega_p[k][3];
}

void update_Tp(Data* data, int k)
{

    double* dTdt = data->dTdt;
    double* Tp = data->Tp;
    double* Sp = data->Sp;
    double** QQ = data->QQ;
    double*** PP = data->PP;
    double** Qr = data->Qr;
    double rho_r = data->rho_r;
    double rho_p = data->rho_p;
    double rho_f = data->rho_f;
    double cr = data->cr;
    double cp = data->cp;
    double cf = data->cf;
    double dt = data->dt;
    double dH = data->dH;

    //compute_Qr(Qr, PP[0][k][2], dH, k);

    dTdt[k] = (23.*QQ[k][2]-16.*QQ[k][1]+5.*QQ[k][0])/(12.*Sp[k]*(rho_r*cr - 1.)) + (23.*Qr[k][2]-16.*Qr[k][1]+5.*Qr[k][0])/(12.*Sp[k]*(rho_p*cp - rho_f*cf));
    Tp[k] += dt*dTdt[k];
    PetscPrintf(PETSC_COMM_WORLD,"Temperature of particle %d: Tp = %f[K] \n", k+1, Tp[k]);
}

void update_Cp(Data* data, int k)
{
    double** dCdt = data->dCdt;
    double** Cp = data->Cp;
    double*** PP = data->PP;
    double dt = data->dt;

    dCdt[k][1] = (23.*(PP[k][0][2]+PP[k][1][2])-16.*(PP[k][0][1]+PP[k][1][1])+5.*(PP[k][0][0]+PP[k][1][0]))/12.;
    Cp[k][0] = 0.;
    Cp[k][1] += dt*dCdt[k][1];
}
*/

double* make1DDoubleArray(int arraySize) {
    double* theArray = (double*) calloc((size_t) arraySize, sizeof(double));
    return theArray;
}

double** make2DDoubleArray(int arraySizeX, int arraySizeY) {
    double** theArray;
    theArray = (double**) malloc(arraySizeX*sizeof(double*));
    for (int ix = 0; ix < arraySizeX; ix++){
        theArray[ix] =(double*) calloc((size_t) arraySizeY, sizeof(double));
    }
    return theArray;
}

double*** make3DDoubleArray(int arraySizeZ, int arraySizeX, int arraySizeY) {
    double*** theArray;
    theArray = (double***) malloc(arraySizeZ*sizeof(double**));
    for (int k = 0; k < arraySizeZ; k++) {
        theArray[k] = (double**) malloc(arraySizeX*sizeof(double*));
    }
    for (int k=0; k < arraySizeZ; k++) {
        for (int ix=0; ix < arraySizeX; ix++) {
            theArray[k][ix] = calloc((size_t) arraySizeY, sizeof(double));
        }
    }
    return theArray;
}

int** make2DIntArray(int arraySizeX, int arraySizeY) {
    int** theArray;
    theArray = (int**) malloc(arraySizeX*sizeof(int*));
    for (int ix = 0; ix < arraySizeX; ix++) {
        theArray[ix] =(int*) calloc( (size_t) arraySizeY, sizeof(int));
    }
    return theArray;
}

int*** make3DIntArray(int numberOfparticles, int arraySizeX, int arraySizeY) {
    int*** theArray;
    theArray = (int***) malloc(numberOfparticles*sizeof(int**));
    for (int k = 0; k < numberOfparticles; k++) {
        theArray[k] = (int**) malloc(arraySizeX*sizeof(int*));
    }
    for (int k=0; k < numberOfparticles; k++) {
        for (int ix=0; ix < arraySizeX; ix++) {
            theArray[k][ix] = calloc( (size_t) arraySizeY, sizeof(int));
        }
    }
    return theArray;
}

void free2Darray(double** array, int dim1){
    for(int i=0; i<dim1; i++){
        free(array[i]);
    }
    free(array);
}

void free3Darray(double*** array, int dim1, int dim2) {
    for (int i = 0; i < dim1; i++) {
        for (int j = 0; j < dim2; j++) {
            free(array[i][j]);
        }
    }
    for (int i = 0; i < dim1; i++) {
        free(array[i]);
    }
    free(array);
}
