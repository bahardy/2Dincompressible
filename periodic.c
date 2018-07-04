//
// Created by Baptiste Hardy on 12/04/18.
//

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
#include "fields_creation.h"
#include "forces.h"
#include "collision.h"
#include "particle_motion.h"


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

    /**------------------------------- Matrix Fields Creation  ------------------------------- **/

    allocate_fields(&data);

    /** ------------------------------- Fields Initialization ------------------------------- **/

    /* Particles position */
    data.xg[0][0] = 15;
    data.xg[0][1] = data.xg[0][0];

    //data.xg[1] = 4.5;
    data.yg[0][0] = data.H;
    data.yg[0][1] = data.yg[0][0];

    data.theta[0][0] = 0;
    data.theta[0][1] = data.theta[0][0];

    data.dp[0] = data.Dp;
    //data.dp[1] = data.Dp;
    data.rp[0] = .5*data.Dp;
    //data.rp[1] = .5*data.Dp;

    data.Up[0][1] = 0;//-data.u_m;
    data.Up[0][0] = data.Up[0][1];
    data.Up[0][2] = data.Up[0][1];

    data.Vp[0][1] = 0.;
    data.Vp[0][0] = data.Vp[0][1];
    data.Vp[0][2] = data.Vp[0][1];


#ifdef TEMP
    /*Initialization of particles temperatures */
    for(int k=0; k<data.Np; k++){
        data.Tp[k] = data.Tp0;
    }

    /* We feed reactants at the inlet */
    data.C0[0] = data.CA0;
    data.C0[1] = data.CB0;
#endif

    for(int k=0; k<data.Np; k++){
#ifdef DISK
        data.Sp[k]=M_PI*data.rp[k]*data.rp[k];
        data.J[k] =(M_PI/2.)*pow(data.rp[k],4);
#endif
#ifdef ELLIPSE
        data.a = 2*Dp;
        data.b = Dp;
        data.Sp[k]=M_PI*data.a*data.b;
        data.J[k] =(M_PI/4.)*data.a*data.b*(pow(data.a,2) + pow(data.b,2));
#endif

    }

    /** -------- Some files creation and data writing-------- **/

    FILE* fichier_data = fopen("results/data.txt", "w");
    writeData(fichier_data, data);
    fclose(fichier_data);
    FILE** fichier_particles = malloc(sizeof(FILE*)*data.Np);
    FILE** fichier_forces = malloc(sizeof(FILE*)*data.Np);
    FILE** fichier_fluxes = malloc(sizeof(FILE*)*data.Np);

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

    FILE* fichier_stat = fopen("results/stats.txt", "w+");
    FILE* fichier_forces_NOCA = fopen("results/forces_NOCA.txt", "w+");

    int m = data.m;
    int n = data.n;
    int Np = data.Np;

    /*Initialization of the mask */
    get_masks(&data);
    get_Us_Vs(&data);

#ifdef WRITE
    /*INITIAL SOLUTION (t=0) AFTER RAMPING */
    if(rank==0){
        writeFields_periodic(&data, 0);
    }
#endif


    /** -------------------------------TIME STEPPING ------------------------------- **/
    data.iter = 1;
    double c1, c2, c3, c4, c5, c6;
    double t = 0;
    double surf = 0;
    double delta;
    double tol = 1e-5;
    double relax;
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
        int k;
        int it = 0;

        /** Check for collisions **/
        collision(&data);

        delta = INFINITY;
        /** Check for collisions **/
        while(delta > tol &&  it < it_max) {
            //collision(&data);
            for (k = 0; k < data.Np; k++) {
                /* Integrate penalization term */
                integrate_penalization_periodic(&data, &surf, k);
#ifdef  MOVE
                /* Velocity - Forces */
                if (t >= data.t_move) {
#ifdef TWO_WAY
                    update_Up(&data, k);
                    data.Up[k][2] = relax*data.Up[k][2] + (1-relax)*Up_old[k];
                    data.Vp[k][2] = relax*data.Vp[k][2] + (1-relax)*Vp_old[k];
                    data.Omega_p[k][2] = relax*data.Omega_p[k][2] + (1-relax)*Omega_p_old[k];
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
            printf("xp = %f \n", Xp_old[0]);
            get_masks(&data);
            get_Us_Vs(&data);
#endif

#ifdef TEMP
            get_Ts(&data);
            get_Cs(&data);
#endif
            get_Ustar_Vstar(&data, data.ramp);

            c1 = fabs(data.xg[0][1] - Xp_old[0]);
            c2 = fabs(data.yg[0][1] - Yp_old[0]);
            c3 = fabs(data.theta[0][1] - theta_old[0])/2*M_PI;
            c4 = fabs(data.Up[0][2] - Up_old[0]);
            c5 = fabs(data.Vp[0][2] - Vp_old[0]);
            c6 = fabs(data.Omega_p[0][2] - Omega_p_old[0]);

            delta = fmax(c1, fmax(c2, fmax(c3, fmax(c4, fmax(c5, c6)))));

            for(k = 0; k<Np; k++)
            {
                Xp_old[k] = data.xg[k][1];
                Yp_old[k] = data.yg[k][1];
                theta_old[k] = data.theta[k][1];

                Up_old[k] = data.Up[k][2];
                Vp_old[k] = data.Vp[k][2];
                Omega_p_old[k] = data.Omega_p[k][2];
            }

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

        update_flow(&data);
        get_ghosts(&data, data.Tm0, data.C0);
        get_vorticity(&data);
        get_tau_periodic(&data);

        diagnostic(&data);

        /** Update quantities **/
        for (k=0; k<Np; k++) {

            data.Fx_coll[k][0] =  data.Fx_coll[k][1];
            data.Fx_coll[k][1] =  data.Fx_coll[k][2];
            data.Fx_coll[k][2] = 0;

            data.Fy_coll[k][0] =  data.Fy_coll[k][1];
            data.Fy_coll[k][1] =  data.Fy_coll[k][2];
            data.Fy_coll[k][2] = 0;

            /* Force along x-direction */
            data.F[k][0] = data.F[k][1]; /* n-2*/
            data.F[k][1] = data.F[k][2]; /* n-1*/

            /* Force along y-direction */
            data.G[k][0] = data.G[k][1];
            data.G[k][1] = data.G[k][2];

            /* Moment along z-direction */
            data.Mz[k][0] = data.Mz[k][1];
            data.Mz[k][1] = data.Mz[k][2];

            data.xg[k][0] = data.xg[k][1];
            data.yg[k][0] = data.yg[k][1];
            data.theta[k][0] = data.theta[k][1];

#ifdef TWO_WAY
            data.Up[k][0] = data.Up[k][1];
            data.Up[k][1] = data.Up[k][2];

            data.Vp[k][0] = data.Vp[k][1];
            data.Up[k][1] = data.Up[k][2];

            data.Omega_p[k][0] = data.Omega_p[k][1];
            data.Omega_p[k][1] = data.Omega_p[k][2];

#endif
        }


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
    free(fichier_forces);
    free(fichier_fluxes);
    free(fichier_particles);

    fclose(fichier_stat);
    fclose(fichier_forces_NOCA);


    /* Free memory */
    free_fields(&data);
    free(Up_old), free(Vp_old), free(Omega_p_old);
    free(Xp_old), free(Yp_old), free(theta_old);


    PetscFinalize();
    return 0;
}


void set_up(Data* data, int argc, char *argv[], int rank)
{
    /* DIMENSIONS */
    data->Dp = 1;
    data->d = 4.;
    data->H = 0.5*data->d;
    data->L = 16.;
    data->h = data->Dp/30;
    data->eps = 0;
#ifdef SMOOTHING
    data->eps = 2*data->h;
#endif
    /* NON-DIMENSIONAL NUMBERS */
    data->Pr = 0.7;
    data->Le = 1; /* Lewis number, ratio between Sc and Prandtl */
    data->Sc = data->Le*data->Pr;
    data->Rep = 40.;
    data->Fr = sqrt(1e3);

    /* FLOW */
    data->u_m = 1.;
    data->g = 0;
#ifdef GRAVITY
    data->g = 1;//9.81;
#endif
    /* ENERGY */
    data->alpha_f = data->nu/data->Pr;
    data->Tm0 = 1; // cup-mixing temperature at the inlet
    data->Tp0 = 0.5; // initial particle temperature

    /* PHYSICAL PARAMETERS */
    data->rho_f = 1.;
    data->rho_p = 10;
    data->rho_r = data->rho_p/data->rho_f;
    data->cp = 1000.;
    data->cf = 1000.;
    data->cr = data->cp/data->cf;

#ifdef SEDIMENTATION
    data->Ga = 1e2;
    data->Rep = data->Ga;
#endif
    data->nu = data->u_m*data->Dp/data->Rep;

    /* SPECIES */
    data->Ns = 2;
    data->Np = 1;
    data->Df = make1DDoubleArray(data->Ns);
    data->Df[0] = data->nu/data->Sc;
    data->Df[1] = data->nu/data->Sc;
    data->dH = 0;
    data->CA0 = 1.;
    data->CB0 = 0.;

    /* GRID */
    data->N = (int) (data->d/data->h);
    data->M = (int) (data->L/data->h);
    data->n = data->N + 2; /*for ghost points */
    data->m = data->M; /*for ghost points */


    /* TIME INTEGRATION */
    data->CFL = 0.1;  /*Courant-Freidrichs-Lewy condition on convective term */
    data->r = .25; /* Fourier condition on diffusive term */
    double dt_CFL = data->CFL*data->h/data->u_m;
    double dt_diff = data->r*data->h*data->h/data->nu;

#ifdef EXPLICIT
    data->ratio_dtau_dt = 1;
#endif
#ifndef EXPLICIT
    data->ratio_dtau_dt = 1e-3;
#endif

    data->dt = fmin(dt_CFL, dt_diff);
    data->dt = 0.005;
    data->dtau = data->ratio_dtau_dt*data->dt;

    data->SORitermax = 100000;
    data->alpha = 1.98;
    data->SORtol = 1e-6;

    if(rank == 0){
        printf("Rep = %f\n", data->Rep);
        printf("ratio L/d = %f \n", data->L/data->d);
        printf("Um = %f\n", data->u_m);
        printf("dt_CFL = %f\n", dt_CFL);
        printf("dt_diff = %f\n", dt_diff);
        printf("CFL condition set to %f h/Umax \n", data->CFL);
        printf("dt = %f\n", data->dt);
        printf("dtau = %f\n", data->dtau);
        printf("ratio dtau/dt = %f \n", data->ratio_dtau_dt);
    }

    /* Writing */

    if (argc >= 3) {
        sscanf(argv[1], "%d", &(data->T_write));
        sscanf(argv[2], "%d", &(data->N_write));
    }
    else{
        data->T_write = 1; /* number of time steps between two writings */
        data->N_write = 50; /* number of times we write in files */
    }

    double Tf = data->N_write*data->T_write*data->dt;
    data->Tf = Tf;
    data->t_move = 0; //data->Tf/10.;
    data->t_transfer = 0;
    data->nKmax = 2;
    data->Kmax = 50; /* number of ramping steps */

    data->ramp = 1;

    if(rank == 0){
        printf("Write every %d * dt \n", data->T_write);
        printf("Write %d times \n", data->N_write);
        printf("Final time : %f \n \n", data->Tf);
    }

}

void compute_Qr(double** Qr, double rate, double dH, int k){

    Qr[k][0] = Qr[k][1];
    Qr[k][1] = Qr[k][2];
    Qr[k][2] = rate*(-dH);
}

void diagnostic(Data* data){
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


void get_ghosts(Data* data, double T0, double* C0)
{
    double** u_n = data->u_n;
    double** v_n = data->v_n;
    double** T_n = data->T_n;
    double*** C = data->C_n;

    int m = data->m;
    int n = data->n;
    int Ns = data->Ns;

    /*Ghost points to impose BC's */
    /* Along y = -H and y = H */
    for (int i=0; i<m; i++){
        /* On u_n */
        /* Bottom and top Wall : slip : du/dn= 0 */
        u_n[i][0] = u_n[i][1];
        u_n[i][n-1] = u_n[i][n-2];

#ifdef TEMP
        /* Walls : adiabatic: dTdn = 0, no mass flux */
        T_n[i][0] = T_n[i][1];
        T_n[i][n-1] = T_n[i][n-2];

        for (int s=0; s<Ns; s++){
            C[s][i][0] = C[s][i][1];
            C[s][i][n-1] = C[s][i][n-2];
        }
#endif
    }

}

void get_masks(Data* data)
{
    double*** chi_S = data->chi_S;
    double*** chi_U = data->chi_U;
    double*** chi_V = data->chi_V;

    double** I_S = data->I_S;
    double** I_U = data->I_U;
    double** I_V = data->I_V;

    double*** Ip_S = data->Ip_S;
    double*** Ip_U = data->Ip_U;
    double*** Ip_V = data->Ip_V;
    double** coloring = data->coloring;
    double* rp = data->rp;
    double dist;
    int m = data->m;
    int n = data->n;
    double h = data->h;
    double L = data->L;
    int Np = data->Np;
    int b1, b2, b3;
    double d1, d2, d3;
    double xU, xV, yU, yV, yS, xS;
    double** xg = data->xg;
    double** yg = data->yg;
    double** theta = data->theta;

    double xloc, yloc, delta;


    for(int i=0; i<m; i++){
        xU = (i+1)*h;
        xV = (i+0.5)*h;
        for(int j=0; j<n; j++){
            yU = (j-0.5)*h;
            yV = j*h;
            yS = yU;
            xS = xV;

            I_S[i][j] = 0; /*Reset the masks */
            I_U[i][j] = 0;
            I_V[i][j] = 0;

            coloring[i][j] = 0;

            /*Go over all the particles */
            for(int k=0; k<Np; k++){
                double xG = fmod(xg[k][1], L);
                double yG = yg[k][1];

#ifdef ELLIPSE
                double x;
		double y;
		/*ELLIPSE*/
		x = xU-Xp_k[k];
		y = yU-yg[k];
                double EU = pow(b,2)*(pow(y*cos(theta[k]),2)+pow(x*sin(theta[k]),2)-y*x*sin(2*theta[k]))
                + pow(a,2)*(pow(x*cos(theta[k]),2)+pow(y*sin(theta[k]),2)+y*x*sin(2*theta[k]))
                - pow(a,2)*pow(b,2);

		x = xV-Xp_k[k];
		y = yV-yG;
                double EV = pow(b,2)*(pow(y*cos(theta[k]),2)+pow(x*sin(theta[k]),2)-y*x*sin(2*theta[k]))
                + pow(a,2)*(pow(x*cos(theta[k]),2)+pow(y*sin(theta[k]),2)+y*x*sin(2*theta[k]))
                - pow(a,2)*pow(b,2);

		x = xS-Xp_k[k];
		y = yS-yG;
                double ES = pow(b,2)*(pow(y*cos(theta[k]),2)+pow(x*sin(theta[k]),2)-y*x*sin(2*theta[k]))
                + pow(a,2)*(pow(x*cos(theta[k]),2)+pow(y*sin(theta[k]),2)+y*x*sin(2*theta[k]))
                - pow(a,2)*pow(b,2);

                Ip_S[k][i][j]=(ES <= 0);
                Ip_U[k][i][j]=(EU <= 0);
                Ip_V[k][i][j]=(EV <= 0);

		coloring[i][j] +=Ip_S[k][i][j];
#endif

#ifdef DISK
                b1 = ((xS-xG)*(xS-xG)+(yS-yG)*(yS-yG)<= (rp[k]+data->eps)*(rp[k]+data->eps));
                b2 = ((xS-(xG+L))*(xS-(xG+L))+(yS-yG)*(yS-yG)<= (rp[k]+data->eps)*(rp[k]+data->eps));
                b3 = ((xS-(xG-L))*(xS-(xG-L))+(yS-yG)*(yS-yG)<= (rp[k]+data->eps)*(rp[k]+data->eps));
                chi_S[k][i][j]= (b1 || b2 || b3);

                b1 = ((xU-xG)*(xU-xG)+(yU-yG)*(yU-yG)<= (rp[k]+data->eps)*(rp[k]+data->eps));
                b2 = ((xU-(xG+L))*(xU-(xG+L))+(yU-yG)*(yU-yG)<= (rp[k]+data->eps)*(rp[k]+data->eps));
                b3 = ((xU-(xG-L))*(xU-(xG-L))+(yU-yG)*(yU-yG)<= (rp[k]+data->eps)*(rp[k]+data->eps));

                chi_U[k][i][j]= (b1 || b2 || b3);

                b1 = ((xV-xG)*(xV-xG)+(yV-yG)*(yV-yG)<= (rp[k]+data->eps)*(rp[k]+data->eps));
                b2 = ((xV-(xG+L))*(xV-(xG+L))+(yV-yG)*(yV-yG)<= (rp[k]+data->eps)*(rp[k]+data->eps));
                b3 = ((xV-(xG-L))*(xV-(xG-L))+(yV-yG)*(yV-yG)<= (rp[k]+data->eps)*(rp[k]+data->eps));

                chi_V[k][i][j]= (b1 || b2 || b3);

#ifndef SMOOTHING
                b1 = ((xS-xG)*(xS-xG)+(yS-yG)*(yS-yG)<= rp[k]*rp[k]);
                b2 = ((xS-(xG+L))*(xS-(xG+L))+(yS-yG)*(yS-yG)<= rp[k]*rp[k]);
                Ip_S[k][i][j]= b1 || b2;

                b1 = ((xU-xG)*(xU-xG)+(yU-yG)*(yU-yG)<= rp[k]*rp[k]);
                b2 = ((xU-(xG+L))*(xU-(xG+L))+(yU-yG)*(yU-yG)<= rp[k]*rp[k]);
                Ip_U[k][i][j]= b1 || b2;

                b1 = ((xV-xG)*(xV-xG)+(yV-yG)*(yV-yG)<= rp[k]*rp[k]);
                b2 = ((xV-(xG+L))*(xV-(xG+L))+(yV-yG)*(yV-yG)<= rp[k]*rp[k]);
                Ip_V[k][i][j]= b1 || b2;

#endif

#ifdef SMOOTHING

                //Smoothing S
                d1 = sqrt((xS-xG)*(xS-xG)+(yS-yG)*(yS-yG));
                d2 = sqrt((xS-(xG+L))*(xS-(xG+L))+(yS-yG)*(yS-yG));
                d3 = sqrt((xS-(xG-L))*(xS-(xG-L))+(yS-yG)*(yS-yG));
                dist =  rp[k] - fmin(d1, fmin(d2,d3));

                if( dist < - data->eps)
                    Ip_S[k][i][j] = 0;
                else if( fabs(dist) <= data->eps)
                    Ip_S[k][i][j] = .5*(1 + dist/data->eps + (1./M_PI)*sin( M_PI* dist/data->eps) );
                else if( dist > data->eps)
                    Ip_S[k][i][j] = 1;

                //Smoothing U
                d1 = sqrt((xU-xG)*(xU-xG)+(yU-yG)*(yU-yG));
                d2 = sqrt((xU-(xG+L))*(xU-(xG+L))+(yU-yG)*(yU-yG));
                d3 = sqrt((xU-(xG-L))*(xU-(xG-L))+(yU-yG)*(yU-yG));
                dist =  rp[k] - fmin(d1, fmin(d2,d3));

                if( dist < - data->eps)
                    Ip_U[k][i][j] = 0;
                else if( fabs(dist) <=data->eps)
                    Ip_U[k][i][j] = .5*(1 + dist/data->eps + (1./M_PI)*sin( M_PI* dist/data->eps) );
                else if( dist > data->eps)
                    Ip_U[k][i][j] = 1;


                //Smoothing V
                d1 = sqrt((xV-xG)*(xV-xG)+(yV-yG)*(yV-yG));
                d2 = sqrt((xV-(xG+L))*(xV-(xG+L))+(yV-yG)*(yV-yG));
                d3 = sqrt((xV-(xG-L))*(xV-(xG-L))+(yV-yG)*(yV-yG));
                dist =  rp[k] - fmin(d1, fmin(d2,d3));

                if( dist < - data->eps)
                    Ip_V[k][i][j] = 0;
                else if( fabs(dist) <= data->eps)
                    Ip_V[k][i][j] = .5*(1 + dist/data->eps + (1./M_PI)*sin( M_PI* dist/data->eps) );
                else if( dist > data->eps)
                    Ip_V[k][i][j] = 1;

#endif

                xloc = xS-xg[k][1];
                yloc = yS-yG;
                delta = atan2(yloc, xloc);
                coloring[i][j] += Ip_S[k][i][j];

                if((int) floor(((delta-theta[k][1])/(M_PI/2.))) % 2 == 0 ){
                    coloring[i][j] = -coloring[i][j];
                }
#endif
                I_S[i][j] += Ip_S[k][i][j];
                I_U[i][j] += Ip_U[k][i][j];
                I_V[i][j] += Ip_V[k][i][j];
            }

            if(I_S[i][j] > 1 || I_U[i][j] > 1 || I_V[i][j] > 1 ){
                PetscPrintf(PETSC_COMM_WORLD, "Collision de particules \n");
                I_S[i][j] = fmin(I_S[i][j], 1);
                I_U[i][j] = fmin(I_U[i][j], 1);
                I_V[i][j] = fmin(I_V[i][j], 1);
            }
        }
    }
}

void get_Cs(Data* data)
{
    double*** Cs = data-> Cs;
    double** Cp = data->Cp;
    double*** chi_S = data->chi_S;
    int m = data->m;
    int n = data->n;
    int Np = data->Np;

    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            Cs[1][i][j] = 0.;
            for(int k=0; k<Np; k++){
                Cs[1][i][j] += chi_S[k][i][j]*Cp[k][1];
            }
        }
    }

}

void get_Ts(Data* data)
{
    double** Ts = data-> Ts;
    double* Tp = data->Tp;
    double*** chi_S = data->chi_S;
    int m = data->m;
    int n = data->n;
    int Np = data->Np;

    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            Ts[i][j] = 0.;
            for(int k = 0; k<Np; k++){
                Ts[i][j]+= chi_S[k][i][j]*Tp[k];

            }
        }
    }
}

void get_Us_Vs(Data* data){

    double** u_s = data-> u_s;
    double** v_s = data-> v_s;
    double** Up = data->Up;
    double** Vp = data->Vp;
    double** Omega_p = data->Omega_p;
    double** xg = data->xg;
    double** yg = data->yg;

    double*** chi_U = data->chi_U;
    double*** chi_V = data->chi_V;

    int m = data->m;
    int n = data->n;
    double L = data->L;
    int Np = data->Np;
    double h = data->h;

    double xV, yU;
    for(int i=0; i<m; i++){
        xV = (i+0.5)*h;
        for(int j=0; j<n; j++){
            yU = (j-0.5)*h;
            u_s[i][j] = 0.;
            v_s[i][j] = 0.;
            for (int k = 0; k<Np; k++){
                u_s[i][j]+= chi_U[k][i][j]*(Up[k][2] - Omega_p[k][2]*(yU-yg[k][1]) );
                v_s[i][j]+= chi_V[k][i][j]*(Vp[k][2] + Omega_p[k][2]*fmod(xV-xg[k][1],L) );
            }
        }
    }
}

void get_Ustar_Vstar(Data* data, double ramp)
{
    int i,j;
    double dt = data->dt;
    double dtau = data->dtau;
    int m = data->m;
    int n = data->n;
    double h = data->h;
    double nu = data->nu;

    double** I_U = data->I_U;
    double** I_V = data->I_V;

    double** u_n = data->u_n;
    double** u_n_1 = data->u_n_1;
    double** u_star = data->u_star;
    double** u_s = data->u_s;

    double** v_n = data->v_n;
    double** v_n_1 = data->v_n_1;
    double** v_star = data->v_star;
    double** v_s = data->v_s;

    double** P = data->P;

    double H_U, H_U_old;
    double H_V, H_V_old;
    double lapU, lapV;
    double dpdx, dpdy;

    double uL, uR, vT, vB;
    double dudxR, dudxL, dudyT, dudyB, dvdxR, dvdxL, dvdyT, dvdyB;

    /* u_star ADAMS-BASHFORTH 2 */
    for (i=0; i<m; i++){
        for (j=1; j<n-1; j++){

            // CONVECTIVE TERM
            /** time n-1 **/
            uR = .5*(u_n_1[i][j] + u_n_1[(i+1+m)%m][j]);
            uL = .5*(u_n_1[(i-1+m)%m][j] + u_n_1[i][j]);
            dudxR = (u_n_1[(i+1+m)%m][j]- u_n_1[i][j])/h;
            dudxL = (u_n_1[i][j]- u_n_1[(i-1+m)%m][j])/h;

            vT = .5*(v_n_1[i][j] + v_n_1[(i+1+m)%m][j]);
            vB = .5*(v_n_1[i][j-1] + v_n_1[(i+1+m)%m][j-1]);
            dudyT = (u_n_1[i][j+1] - u_n_1[i][j])/h;
            dudyB = (u_n_1[i][j] - u_n_1[i][j-1])/h;

            H_U_old = .5*(uR*dudxR + uL*dudxL) + .5*(vT*dudyT + vB*dudyB);

            /** time n **/

            uR = .5*(u_n[i][j] + u_n[(i+1+m)%m][j]);
            uL = .5*(u_n[(i-1+m)%m][j] + u_n[i][j]);
            dudxR = (u_n[(i+1+m)%m][j]- u_n[i][j])/h;
            dudxL = (u_n[i][j]- u_n[(i-1+m)%m][j])/h;

            vT = .5*(v_n[i][j] + v_n[(i+1+m)%m][j]);
            vB = .5*(v_n[i][j-1] + v_n[(i+1+m)%m][j-1]);
            dudyT = (u_n[i][j+1] - u_n[i][j])/h;
            dudyB = (u_n[i][j] - u_n[i][j-1])/h;

            H_U = .5*(uR*dudxR + uL*dudxL) + .5*(vT*dudyT + vB*dudyB);


            // LAPLACIAN
            lapU = (u_n[(i+1+m)%m][j]+u_n[(i-1+m)%m][j]+u_n[i][j+1]+u_n[i][j-1]-4.*u_n[i][j])/(h*h);

            // PRESSURE TERM
            dpdx = (P[(i+1+m)%m][j]-P[i][j])/h;

#ifdef EXPLICIT
            //EXPLICIT VERSION
            u_star[i][j] = u_n[i][j] + dt*(-1.5*H_U + 0.5*H_U_old - dpdx + nu*lapU - ramp*I_U[i][j]*(u_n[i][j] - u_s[i][j])/dtau);

#else
            //IMPLICIT VERSION
            u_star[i][j] = (u_n[i][j] + dt*(-1.5*H_U + 0.5*H_U_old - dpdx + nu*lapU) + (dt/dtau)*ramp*I_U[i][j]*u_s[i][j])/(1.+ramp*I_U[i][j]*dt/dtau);
#endif

        }
    }

    /* v_star  ADAMS-BASHFORTH 2 */
    for (i=0; i<m; i++){
        for (j=1; j<n-2; j++){

            // CONVECTIVE TERM

            uR = .5*(u_n_1[i][j] + u_n_1[i][j+1]);
            uL = .5*(u_n_1[(i-1+m)%m][j] + u_n_1[(i-1+m)%m][j+1]);
            dvdxR = (v_n_1[(i+1+m)%m][j]- v_n_1[i][j])/h;
            dvdxL = (v_n_1[i][j]- v_n_1[(i-1+m)%m][j])/h;

            vT = .5*(v_n_1[i][j] + v_n_1[i][j+1]);
            vB = .5*(v_n_1[i][j] + v_n_1[i][j-1]);
            dvdyT = (v_n_1[i][j+1] - v_n_1[i][j])/h;
            dvdyB = (v_n_1[i][j] - v_n_1[i][j-1])/h;

            H_V_old = .5*(uR*dvdxR + uL*dvdxL) + .5*(vT*dvdyT + vB*dvdyB);

            uR = .5*(u_n[i][j] + u_n[i][j+1]);
            uL = .5*(u_n[(i-1+m)%m][j] + u_n[(i-1+m)%m][j+1]);
            dvdxR = (v_n[(i+1+m)%m][j]- v_n[i][j])/h;
            dvdxL = (v_n[i][j]- v_n[(i-1+m)%m][j])/h;

            vT = .5*(v_n[i][j] + v_n[i][j+1]);
            vB = .5*(v_n[i][j] + v_n[i][j-1]);
            dvdyT = (v_n[i][j+1] - v_n[i][j])/h;
            dvdyB = (v_n[i][j] - v_n[i][j-1])/h;

            H_V = .5*(uR*dvdxR + uL*dvdxL) + .5*(vT*dvdyT + vB*dvdyB);

            // LAPLACIAN
            lapV = (v_n[(i+1+m)%m][j]+v_n[(i-1+m)%m][j]+v_n[i][j+1]+v_n[i][j-1]-4.*v_n[i][j])/(h*h);

            // PRESSURE TERM
            dpdy = (P[i][j+1]-P[i][j])/h;


#ifdef EXPLICIT
            //EXPLICIT VERSION
            v_star[i][j] = v_n[i][j] + dt*(-1.5*H_V + 0.5*H_V_old - dpdy + nu*lapV - ramp*I_V[i][j]*(v_n[i][j] - v_s[i][j])/dtau);
#else
            //IMPLICIT VERSION
            v_star[i][j] = (v_n[i][j] + dt*(-1.5*H_V + 0.5*H_V_old - dpdy + nu*lapV) + (dt/dtau)*ramp*I_V[i][j]*v_s[i][j])/(1.+ramp*I_V[i][j]*dt/dtau);
#endif
        }
    }

}

void update_flow(Data* data) {

    int m = data->m;
    int n = data->n;
    double dt = data->dt;
    double h = data->h;
    //double H = data->H;

    double **u_new = make2DDoubleArray(m, n);
    double **v_new = make2DDoubleArray(m, n);

    double **u_n_1 = data->u_n_1;
    double **u_n = data->u_n;
    double **v_n_1 = data->v_n_1;
    double **v_n = data->v_n;
    double **P = data->P;

    double **u_star = data->u_star;
    double **v_star = data->v_star;
    double **phi = data->phi;

    int i, j;
    /* Set boundary for u_new, v_new */

    for (i = 0; i < m; i++) {
        for (j = 0; j < n - 1; j++) {
            u_new[i][j] = u_star[i][j] - dt * (phi[(i+1 + m) % m][j] - phi[i][j]) / h;
            v_new[i][j] = v_star[i][j] - dt * (phi[i][j+1] - phi[i][j]) / h;
            P[i][j] += phi[i][j];
        }
    }

#ifdef TEMP
    /* TEMPERATURE AND SPECIES */

    double **T_new = make2DDoubleArray(m, n);
    double **T_n = data->T_n;
    double **T_n_1 = data->T_n_1;
    double **Ts = data->Ts;
    double ***Cs = data->Cs;

    int Ns = data->Ns;
    double ***C_new = make3DDoubleArray(Ns, m, n);
    double ***C_n = data->C_n;
    double ***C_n_1 = data->C_n_1;

    double alpha_f = data->alpha_f;
    double *Df = data->Df;
    double **I_S = data->I_S;
    double ramp = data->ramp;
    double dtau = data->dtau;

    double Uij, Uij_old, Vij, Vij_old;
    double H_T, H_T_old, lapT;
    double H_C, H_C_old, lapC;

    for (i = 1; i < m - 1; i++) {
        for (j = 1; j < n - 1; j++) {

            // Need to have a value for velocities at cell center
            Uij = .5 * (u_n[i][j] + u_n[i - 1][j]); /* average of the two neighbour points along x-direction */
            Vij = .5 * (v_n[i][j] + v_n[i - 1][j]); /* average of the two neighbour points along y-direction */
            Uij_old = .5 * (u_n_1[i][j] + u_n_1[i - 1][j]);
            Vij_old = .5 * (v_n_1[i][j] + v_n_1[i - 1][j]);

            // Advective terms
            H_T_old = Uij_old * (T_n_1[i + 1][j] - T_n_1[i - 1][j]) / (2. * h) +
                      Vij_old * (T_n_1[i][j + 1] - T_n_1[i][j - 1]) / (2. * h);
            H_T = Uij * (T_n[i + 1][j] - T_n[i - 1][j]) / (2. * h) + Vij * (T_n[i][j + 1] - T_n[i][j - 1]) / (2. * h);
            // Laplacian
            lapT = (T_n[i + 1][j] + T_n[i - 1][j] + T_n[i][j + 1] + T_n[i][j - 1] - 4. * T_n[i][j]) / (h * h);

            T_new[i][j] = (T_n[i][j] + dt * (-1.5 * H_T + 0.5 * H_T_old
                                             + alpha_f * lapT
                                             + ramp * I_S[i][j] * Ts[i][j] / dtau)) /
                          (1. + ramp * I_S[i][j] * dt / dtau);

            for (int s = 0; s < Ns; s++) {
                // Advective terms
                H_C_old = Uij_old * (C_n_1[s][i + 1][j] - C_n_1[s][i - 1][j]) / (2. * h) +
                          Vij_old * (C_n_1[s][i][j + 1] - C_n_1[s][i][j - 1]) / (2. * h);
                H_C = Uij * (C_n[s][i + 1][j] - C_n[s][i - 1][j]) / (2. * h) +
                      Vij * (C_n[s][i][j + 1] - C_n[s][i][j - 1]) / (2. * h);
                // Laplacian
                lapC = (C_n[s][i + 1][j] + C_n[s][i - 1][j] + C_n[s][i][j + 1] + C_n[s][i][j - 1] - 4. * C_n[s][i][j]) / (h * h);

                C_new[s][i][j] = (C_n[s][i][j] + dt * (-1.5 * H_C + 0.5 * H_C_old
                                                       + Df[s] * lapC
                                                       + ramp * I_S[i][j] * Cs[s][i][j] / dtau)) /
                                 (1. + ramp * I_S[i][j] * dt / dtau);
            }

        }
    }

    free2Darray(T_n_1, m);
    free3Darray(C_n_1, Ns, m);

    data->T_n_1 = T_n;
    data->T_n = T_new;

    data->C_n_1 = C_n;
    data->C_n = C_new;

#endif

    free2Darray(u_n_1, m);
    free2Darray(v_n_1, m);

    data->u_n_1 = u_n;
    data->u_n = u_new;

    data->v_n_1 = v_n;
    data->v_n = v_new;

}


void get_vorticity(Data* data){
    int m = data->m;
    int n = data->n;
    double h = data->h;

    int i,j;

    double** omega = data->omega;
    double** u_n = data->u_n;
    double** v_n = data->v_n;

    for (i = 0; i < m; i++) {
        for (j = 1; j < n - 2; j++) {
            omega[i][j] = (v_n[(i+1+m)%m][j] - v_n[i][j]) / h - (u_n[i][j + 1] - u_n[i][j]) / h;
        }
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

