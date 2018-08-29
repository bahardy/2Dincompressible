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
#include "ghost_points.h"


int main(int argc, char *argv[]){

    PetscInitialize(&argc, &argv, 0, 0);
    int rank, nbproc;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank );
    MPI_Comm_size(PETSC_COMM_WORLD, &nbproc);

    struct stat st = {0};
    char* folder = "results";
    if (stat(folder, &st) == -1) {
        mkdir(folder, 0700);
    }

    /**------------------------------- DATA BASE CREATION ------------------------------- **/
    Data data;
    set_up(&data, argc, argv, rank);

    FILE* fichier_data = NULL;
    char fileData[30];
    strcpy(fileData, folder);
    strcat(fileData, "/data.txt");
    fichier_data = fopen(fileData, "w");
    writeData(fichier_data, data);
    fclose(fichier_data);

    /** -------- Some files creation and data writing-------- **/

    FILE** fichier_particles = malloc(sizeof(FILE*)*data.Np);
    FILE** fichier_forces = malloc(sizeof(FILE*)*data.Np);
    FILE** fichier_fluxes = malloc(sizeof(FILE*)*data.Np);

    FILE* fichier_stat = NULL;
    char fileStat[30];
    strcpy(fileStat, folder);
    strcat(fileStat, "/stats.txt");
    fichier_stat = fopen(fileStat, "w");

    FILE* fichier_NOCA = NULL;
    char fileNOCA[30];
    strcpy(fileNOCA, folder);
    strcat(fileNOCA, "/forces_NOCA.txt");
    fichier_NOCA = fopen(fileNOCA, "w");

    for(int k = 0; k<data.Np; k++)
    {
        char K[10];
        sprintf(K, "%d", k);

        char fileParticle[30];
        strcpy(fileParticle, folder);
        strcat(fileParticle, "/particle");
        strcat(fileParticle, "-");
        strcat(fileParticle, K);
        strcat(fileParticle, ".txt");
        fichier_particles[k] = fopen(fileParticle, "w+");

        char fileForces[30];
        strcpy(fileForces, folder);
        strcat(fileForces, "/forces");
        strcat(fileForces, "-");
        strcat(fileForces, K);
        strcat(fileForces, ".txt");
        fichier_forces[k] = fopen(fileForces, "w+");

        char fileFluxes[30];
        strcpy(fileFluxes, folder);
        strcat(fileFluxes, "/fluxes");
        strcat(fileFluxes, "-");
        strcat(fileFluxes, K);
        strcat(fileFluxes, ".txt");
        fichier_fluxes[k] = fopen(fileFluxes, "w+");
    }

    /** ------------------------------- FIELDS INITIALIZATION ------------------------------- **/
    allocate_fields(&data);
    initialize_fields(&data);
    get_masks(&data);
    get_normal(&data);

    /** ----- BOUNDARY CONDITION -----------**/
    get_ghosts(&data, data.T0, data.C0);


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
        writeFields(&data, folder, 0);
    }
#endif

    /** -------------------------------TIME STEPPING ------------------------------- **/
    double t = 0;
    double surf = 0.;
    int k = 0;

    data.ramp = 1;
    data.iter = 1;

    while(t < data.Tf){

        PetscPrintf(PETSC_COMM_WORLD, "\n \n BEGIN iter %d : t = %f \n", data.iter, t);
        for (k = 0; k<data.Np; k++) {
            PetscPrintf(PETSC_COMM_WORLD, "Position of the center of mass of particle %d: (x,y) = (%f,%f) \n", k + 1,
                        data.xg[k][2], data.yg[k][2]);
            PetscPrintf(PETSC_COMM_WORLD, "Angle: theta  = %f \n", data.theta[k][2]);
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
#ifdef EE
                if(t > data.t_coupling) {
                    update_Up(&data, k);
                }
                update_Xp(&data, k);
#endif
#ifdef AB3
                update_Xp(&data, k);
                if(t > data.t_coupling) {
                    update_Up(&data, k);
                }

#endif
#ifdef LF
		update_Xp(&data, k);
		if(t > data.t_coupling) {
                    update_Up(&data, k);
                }
#endif
            }
#endif
#ifdef  TEMP
            /*Temperature - Species - Fluxes */
#ifndef INTRAPARTICLE
//            if(t > data.t_transfer)
//            {
//                update_Tp(&data, k);
//                update_Cp(&data, k);
//            }
#endif
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
        get_normal(&data);
        get_Us_Vs(&data);
#endif

#ifdef TEMP
#ifndef INTRAPARTICLE
        get_Ts(&data);
        //get_Cs(&data);
#else
        get_conductivity(&data);
        get_diffusivity(&data);
#endif

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
        get_ghosts(&data, data.T0, data.C0);
        get_vorticity(&data);
        get_tau(&data);

        compute_forces_NOCA(&data, fichier_NOCA, 1, data.m-2, 1, data.n-2);
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
                writeFields(&data, folder, data.iter);
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
    fclose(fichier_NOCA);

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
        data->xg[k][1] = data->xg[k][2];

        data->yg[k][0] = data->yg[k][1];
        data->yg[k][1] = data->yg[k][2];

        data->theta[k][0] = data->theta[k][1];
        data->theta[k][1] = data->theta[k][2];
#ifndef AB3
	data->Up[k][0] = data->Up[k][1];
        data->Up[k][1] = data->Up[k][2];

	data->Vp[k][0] = data->Vp[k][1];
	data->Vp[k][1] = data->Vp[k][2];

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
