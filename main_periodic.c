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
#include "flow_solver_periodic.h"


int main(int argc, char *argv[]){

    PetscInitialize(&argc, &argv, 0, 0);
    int rank, nbproc;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank );
    MPI_Comm_size(PETSC_COMM_WORLD, &nbproc);

    struct stat st = {0};
    if (stat("results", &st) == -1) {
        mkdir("results", 0700);
    }

    /**------------------------------- DATA BASE CREATION ------------------------------- **/
    Data data;
    set_up_periodic(&data, argc, argv, rank);
    FILE* fichier_data = fopen("results/data.txt", "w");
    writeData(fichier_data, data);
    fclose(fichier_data);

    /** -------- Some files creation -------- **/

    FILE** fichier_particles = malloc(sizeof(FILE*)*data.Np);
    FILE** fichier_forces = malloc(sizeof(FILE*)*data.Np);
    FILE** fichier_fluxes = malloc(sizeof(FILE*)*data.Np);
    FILE* fichier_stat = fopen("results/stats.txt", "w+");

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



    /**------------------------------- FIELDS INITIALIZATION  ------------------------------- **/

    allocate_fields(&data);
    initialize_fields_periodic(&data);
    get_masks(&data);
    get_Us_Vs(&data);

#ifdef WRITE
    /*INITIAL SOLUTION (t=0) */
    if(rank==0){
        writeFields_periodic(&data, 0);
    }
#endif


    /** -------------------------------TIME STEPPING ------------------------------- **/
    data.iter = 1;
    double t = 0;
    int Np = data.Np;

    double surf = 0;
    double tol = 1e-5;
    double delta;
    double relax;
    int it;
    int it_max;

#ifndef ITERATIVE
    it_max = 1;
    relax = 1;
#endif
#ifdef ITERATIVE
    it_max = 100;
    relax = 0.5;
#endif


    double* Up_old =  make1DDoubleArray(Np);
    double* Vp_old =  make1DDoubleArray(Np);
    double* Omega_p_old =  make1DDoubleArray(Np);

    double* Xp_old = make1DDoubleArray(Np);
    double* Yp_old = make1DDoubleArray(Np);
    double* theta_old = make1DDoubleArray(Np);

    int i, j;
    for (int K = 0; K< Np; K++) {
        Xp_old[K] = data.xg[K][0];
        Yp_old[K] = data.yg[K][0];
        theta_old[K] = data.theta[K][0];

        Up_old[K] = data.Up[K][0];
        Vp_old[K] = data.Vp[K][0];
        Omega_p_old[K] = data.Omega_p[K][0];
    }


    while(t < data.Tf){

        PetscPrintf(PETSC_COMM_WORLD, "\n \n BEGIN iter %d : t = %f \n", data.iter, t);

        /** --- SOLVE SOLID PHASE --- */
        int k = 0;
        it = 0;
        delta = INFINITY;

        /** Check for collisions **/
        while(delta > tol &&  it < it_max) {
            collision(&data);

            for (k = 0; k < data.Np; k++) {
                integrate_penalization_periodic(&data, &surf, k);
#ifdef  MOVE
                if (t >= data.t_move) {
#ifdef TWO_WAY
                    update_Up(&data, k);
                    relax_Up(&data, relax, Up_old, Vp_old, Omega_p_old, k);
#endif
                    update_Xp(&data, k);
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
            }

            /** --- SOLVE FLUID PHASE --- */
#ifdef  MOVE
            printf("xp = %f \n", data.xg[0][1]);
            get_masks(&data);
            get_Us_Vs(&data);
#endif

#ifdef TEMP
            get_Ts(&data);
            get_Cs(&data);
#endif
            get_Ustar_Vstar(&data, data.ramp);

            delta = check_convergence(&data, Xp_old, Yp_old, theta_old, Up_old, Vp_old, Omega_p_old);

            it++;
        }

        PetscPrintf(PETSC_COMM_WORLD, "Fluid-solid coupling achievd after %d iterations. Delta = %f \n", it, delta);
        PetscPrintf(PETSC_COMM_WORLD,"Position of the center of mass of particle %d: (x,y) = (%f,%f) \n", k+1, data.xg[0][1], data.yg[0][1]);
        PetscPrintf(PETSC_COMM_WORLD,"Angle: theta  = %f \n", data.theta[0][1]);

        compute_forces_fluxes(&data, 0);

        clock_t t_init = clock();
        poisson_solver_periodic(&data, rank, nbproc);
        clock_t t_final = clock();
        double t_Poisson = ((double) (t_final - t_init)) / CLOCKS_PER_SEC;
        PetscPrintf(PETSC_COMM_WORLD, "Poisson solver took %f seconds \n", t_Poisson);
        poisson_residual_periodic(&data);
#ifdef TEMP
        update_scalars(&data);
#endif
        update_flow(&data);
        get_ghosts(&data, data.Tm0, data.C0);
        get_vorticity(&data);
        get_tau_periodic(&data);

        diagnostic(&data);
        update_quantities(&data);

#ifdef WRITE
        if(rank == 0){
            fprintf(fichier_stat, "%3.13e \t  %3.13e \t  %3.13e \n", data.CFL_max, data.Reh_max, data.Reh_omega_max);
            fflush(fichier_stat);

            for (k = 0; k< data.Np; k++) {
                writeForces(&data, fichier_forces, k);
                writeParticle(&data, fichier_particles, k);
                writeFluxes(&data, fichier_fluxes, k);
            }

            if(data.iter % data.T_write == 0){
                writeFields_periodic(&data, data.iter);
            }
        }
#endif
        /* Increment time */
        t += data.dt;
        data.iter ++;

    }


    for (int k = 0; k< data.Np; k++)
    {
        fclose(fichier_forces[k]);
        fclose(fichier_fluxes[k]);
        fclose(fichier_particles[k]);
    }
    fclose(fichier_stat);

    free(fichier_forces);
    free(fichier_fluxes);
    free(fichier_particles);
    
    /* Free memory */
    free_fields(&data);
    free(Up_old), free(Vp_old), free(Omega_p_old);
    free(Xp_old), free(Yp_old), free(theta_old);


    PetscFinalize();
    return 0;
}

double check_convergence(Data* data, double* Xp_old, double* Yp_old, double* theta_old, double* Up_old, double* Vp_old, double* Omega_p_old)
{
    double c1, c2, c3, c4, c5, c6;
    c1 = fabs(data->xg[0][1] - Xp_old[0]);
    c2 = fabs(data->yg[0][1] - Yp_old[0]);
    c3 = fabs(data->theta[0][1] - theta_old[0])/2*M_PI;
    c4 = fabs(data->Up[0][2] - Up_old[0]);
    c5 = fabs(data->Vp[0][2] - Vp_old[0]);
    c6 = fabs(data->Omega_p[0][2] - Omega_p_old[0]);

    double delta = fmax(c1, fmax(c2, fmax(c3, fmax(c4, fmax(c5, c6)))));

    for(int k = 0; k<data->Np; k++)
    {
        Xp_old[k] = data->xg[k][1];
        Yp_old[k] = data->yg[k][1];
        theta_old[k] = data->theta[k][1];

        Up_old[k] = data->Up[k][2];
        Vp_old[k] = data->Vp[k][2];
        Omega_p_old[k] = data->Omega_p[k][2];
    }

    return delta;
}

void diagnostic(Data* data)
{
    double** u_n = data->u_n;
    double** v_n = data->v_n;
    double** omega = data->omega;
    double** Reh = data->Reh;
    double** Reh_omega = data->Reh_omega;
    double** CFL_array = data->CFL_array;
    double dt = data->dt;
    double h = data->h;
    double nu = data->nu;
    double m = data->m;
    double n = data->n;

    data->Reh_max = 0.;
    data->Reh_omega_max = 0.;
    data->CFL_max =0.;

    int i, j;
    for (i=0; i<m; i++) {
        for (j = 1; j<n-1; j++) {
            CFL_array[i][j] = (fabs(u_n[i][j]) + fabs(v_n[i][j]))*dt/h;
            Reh[i][j] = (fabs(u_n[i][j]) + fabs(v_n[i][j]))*h/nu;
            Reh_omega[i][j] = fabs(omega[i][j])*h*h/n;
            if(CFL_array[i][j] > data->CFL_max){
                data->CFL_max = CFL_array[i][j];
            }
            if(Reh[i][j] > data->Reh_max){
                data->Reh_max = Reh[i][j];
            }
            if(Reh_omega[i][j] > data->Reh_omega_max){
                data->Reh_omega_max = Reh_omega[i][j];
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

double*** make3DDoubleArray(int numberOfparticles, int arraySizeX, int arraySizeY) {
    double*** theArray;
    theArray = (double***) malloc(numberOfparticles*sizeof(double**));
    for (int k = 0; k < numberOfparticles; k++) {
        theArray[k] = (double**) malloc(arraySizeX*sizeof(double*));
    }
    for (int k=0; k < numberOfparticles; k++) {
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
        free(array[i]);
    }
    free(array);
}

