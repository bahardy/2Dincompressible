#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <petsc.h>
#include <sys/stat.h>

#include "main.h"
#include "write.h"
#include "fields_creation.h"
#include "poisson.h"
#include "forces.h"
#include "collision.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

//#define RAMPING
//#define EXPLICIT
#define WRITE
#define DISK
#define SLIP
//#define GRAVITY
#define SMOOTHING


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
    double t = 0;
    data.ramp = 1;
    data.iter = 0;
    double surf = 0.;


    /**------------------------------- Matrix Fields Creation  ------------------------------- **/
    initialize_fields(&data);

    /** ------------------------------- Fields Initialization ------------------------------- **/

    /* Particles position */
    data.xg[0] = 3;
    data.yg[0] = data.H;
    data.dp[0] = data.Dp;
    data.rp[0] = .5*data.Dp;
    data.theta[0] = 0; // M_PI/10.

    for(int k=0; k<data.Np; k++){
#ifdef DISK
        data.Sp[k]=M_PI*data.rp[k]*data.rp[k];
        //data.II[k]=(data.dp[k]*data.dp[k])/8.; /* IN 2-D !! */
        data.J[k] =(M_PI/2.)*pow(data.rp[k],4);
#endif
#ifdef ELLIPSE
        data.a = 2*Dp;
        data.b = Dp;
        data.Sp[k]=M_PI*data.a*data.b;
        data.II[k]=.25*(pow(data.a,2) + pow(data.b,2));
        data.J[k] =(M_PI/4.)*data.a*data.b*(pow(data.a,2) + pow(data.b,2));
#endif
    }

    /** -------- Some files creation and data writing-------- **/

    FILE* fichier_data = fopen("results/data.txt", "w");
    writeData(fichier_data, data);
    fclose(fichier_data);

    FILE* fichier_stat = fopen("results/stats.txt", "w");
    FILE* fichier_forces = fopen("results/forces.txt", "w");
    FILE* fichier_fluxes = fopen("results/fluxes.txt", "w");
    FILE* fichier_particle = fopen("results/particle.txt", "w");

    /** ------------------------------- INITIALIZATION of the domain ------------------------------- **/

    /* VELOCITY : horizontal flow Um  */
    for(int i=0; i<data.m; i++){
        for(int j=0; j<data.n; j++){
            data.u_n[i][j] = 0;//data.u_m;
            data.u_n_1[i][j] = data.u_n[i][j];
            data.u_star[i][j] = data.u_n[i][j];
            data.T_n[i][j] = data.Tm0;
            data.T_n_1[i][j] = data.T_n[i][j];
            /* v_n is initially at zero */
        }
    }

    /*Initialization of particles temperatures */
    for(int k=0; k<data.Np; k++){
        data.Tp[k] = data.Tp0;
    }

    for(int k=0; k<data.Np; k++){
        data.Up[k][2] = 0.5*data.u_m;
        data.Up[k][1] = data.Up[k][2];
        data.Up[k][0] = data.Up[k][2];

        data.Vp[k][2] = data.u_m;
        data.Vp[k][1] = data.Vp[k][2];
        data.Vp[k][0] = data.Vp[k][2];

    }

    /** ----- BOUNDARY CONDITION -----------**/

    get_ghosts(&data, data.Tm0, data.C0);

    /*Initialization of the mask */
    get_masks(&data);

    /** ------------------------------- RAMPING ------------------------------- **/

    /* Starting from a flow without the particles (i.e. ramp = 0 everywhere) and
     then progressively increasing ramp until ramp*chi = 1 at the position of the ellipse.
     At the end of the loop we obtain the steady solution of the flow for a given
     fixed position of the ellipse. */
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

    data.iter = 1;

    // we feed reactants at the inlet
    data.C0[0] = data.CA0;
    //data.C0[1] = data.CB0;

    while(t < data.Tf){

        PetscPrintf(PETSC_COMM_WORLD, "\n \n BEGIN iter %d : t = %f \n", data.iter, t);

        /** --- SOLVE SOLID PHASE --- */
        int flag_out = 0;
        int k;

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
        /* Deduce solid velocity field */
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

        update_flow(&data);

        get_ghosts(&data, data.Tm0, data.C0);
        get_vorticity(&data);

        diagnostic(&data);


#ifdef WRITE
        if(rank == 0){
            writeStatistics(&data, fichier_stat);
            writeForces(&data, fichier_forces);
            writeParticle(&data, fichier_particle);
            writeFluxes(&data, fichier_fluxes);

            if(data.iter % data.T_write == 0){
                writeFields(&data, data.iter);
            }
        }
#endif
        /* Increment time */
        t += data.dt;
        data.iter ++;

    }
    fclose(fichier_stat);
    fclose(fichier_forces);
    fclose(fichier_fluxes);
    fclose(fichier_particle);

    free_fields(&data);

    PetscFinalize();
    return 0;
}

void set_up(Data* data, int argc, char* argv[], int rank)
{
    /* DIMENSIONS */
    data->Dp = 1.;
    data->d = 4.;
    data->H = 0.5*data->d;
    data->L = 8.;
    data->h = data->Dp/30;
    data->eps = 0;
#ifdef SMOOTHING
    data->eps = 2.*data->h;
#endif
    /* NON-DIMENSIONAL NUMBERS */
    data->Pr = 0.7;
    data->Le = 1; /* Lewis number, ratio between Sc and Prandtl */
    data->Sc = data->Le*data->Pr;
    data->Rep = 40.;
    data->Fr = sqrt(1e3);

    /* FLOW */
    data->u_m = 1.;
    //data->u_max = 1.5*data->u_m;
    //data->nu = 0.03926;
    data->nu = data->u_m*data->Dp/data->Rep;
    data->g = 0;
#ifdef GRAVITY
    data->g = 9.81;
#endif
    /* ENERGY */
    data->alpha_f = data->nu/data->Pr;
    data->Tm0 = 1; // cup-mixing temperature at the inlet
    data->Tp0 = 0.5; // initial particle temperature

    /* PHYSICAL PARAMETERS */
    data->rho_f = 1.;
    data->rho_p = 100.;
    data->rho_r = data->rho_p/data->rho_f;
    data->cp = 1000.;
    data->cf = 1000.;
    data->cr = data->cp/data->cf;

    /* SPECIES */
    data->Ns = 1;
    data->Np = 1;
    data->Df = make1DDoubleArray(data->Ns);
    for (int i = 0; i < data->Ns; i++)
    {
        data->Df[i] = data->nu/data->Sc;
    }
//    data->Df[0] = data->nu/data->Sc;
//    data->Df[1] = data->nu/data->Sc;
    data->dH = 0;
    data->CA0 = 1.;
    data->CB0 = 0.;

    /* GRID */
    data->N = (int) (data->d/data->h);
    data->M = (int) (data->L/data->h);
    data->n = data->N + 2; /*for ghost points */
    data->m = data->M + 2; /*for ghost points */


    /* TIME INTEGRATION */
    data->CFL = 0.1; /*Courant-Freidrichs-Lewy condition on convective term */
    data->r = .25; /* Fourier condition on diffusive term */
    double dt_CFL = data->CFL*data->h/data->u_m;
    double dt_diff = data->r*data->h*data->h/data->nu;

    data->ratio_dtau_dt = 1e-3;
    data->dt = fmin(dt_CFL, dt_diff);
    data->dtau = data->ratio_dtau_dt*data->dt;

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
        printf("smoothing : eps/h = %f \n", data->eps/data->h);
    }

    /* Writing */

    if (argc >= 3) {
        sscanf(argv[1], "%d", &(data->T_write));
        sscanf(argv[2], "%d", &(data->N_write));
    }
    else{
        data->T_write = 1; /* number of time steps between two writings */
        data->N_write = 3; /* number of times we write in files */
    }

    double Tf = data->N_write*data->T_write*data->dt;
    data->Tf = Tf;
    data->t_move = 0; //data->Tf/10.;
    data->t_transfer = 3.;
    data->nKmax = 2;
    data->Kmax = 50; /* number of ramping steps */


    if(rank == 0){
        printf("Write every %d * dt \n", data->T_write);
        printf("Write %d times \n", data->N_write);
        printf("Final time : %f \n \n", data->Tf);
    }
}

void compute_Qr(double** Qr, double rate, double dH, int k)
{

    Qr[k][0] = Qr[k][1];
    Qr[k][1] = Qr[k][2];
    Qr[k][2] = rate*(-dH);
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
#ifdef SLIP
        /* Bottom and top Wall : slip : du/dn= 0 */
        u_n[i][0] = u_n[i][1];
        u_n[i][n-1] = u_n[i][n-2];
#endif
#ifndef SLIP
        /* No-slip (channel) */
        u_n[i][0] = -0.2*(u_n[i][3] - 5.*u_n[i][2] + 15.*u_n[i][1]);
        u_n[i][n-1] = -0.2*(u_n[i][n-4] - 5.*u_n[i][n-3] + 15.*u_n[i][n-2]);
#endif

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

    /* Along x = 0 and x = L */
    for (int j = 0; j<n; j++){
        /* Inflow : horizontal flow --> v_n = 0 */
        v_n[0][j] = -0.2*(v_n[3][j] - 5.*v_n[2][j] + 15.*v_n[1][j]);
        /* Natural outflow : dV/dx =0 */
        v_n[m-1][j] = v_n[m-2][j];

#ifdef TEMP
        /* On T_n and C */
        /* Inflow : T_n uniform  */
        T_n[0][j] = -0.2*(T_n[3][j]-5.*T_n[2][j]+15.*T_n[1][j]-16.*T0);
        T_n[m-1][j] = (7.*T_n[m-2][j]-5.*T_n[m-3][j]+T_n[m-4][j])/3.;

//        /* Inflow : CA = CA0; CB = CB0 */
//        C[0][0][j] = -0.2*(C[0][3][j]-5.*C[0][2][j]+15.*C[0][1][j]-16.*C0[0]);
//        C[1][0][j] = -0.2*(C[1][3][j]-5.*C[1][2][j]+15.*C[1][1][j]-16.*C0[1]);

        /*Outflow : We cancel axial dispersion d2T/dx2 = 0; d2C/dx2 = 0; */
        for(int s=0; s<Ns; s++){
            C[s][0][j] = -0.2*(C[s][3][j]-5.*C[s][2][j]+15.*C[s][1][j]-16.*C0[s]);
            C[s][m-1][j] = (7.*C[s][m-2][j]-5.*C[s][m-3][j]+C[s][m-4][j])/3.;
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
    double* xg = data->xg;
    double* yg = data->yg;
    double* theta = data->theta;
    double* rp = data->rp;
    double d;
    int m = data->m;
    int n = data->n;
    double h = data->h;
    int Np = data->Np;

    double xU, xV, yU, yV, yS, xS;
    double xloc, yloc, delta;

    for(int i=0; i<m; i++){
        xU = i*h;
        xV = (i-0.5)*h;
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

#ifdef ELLIPSE
                double x;
		double y;
		/*ELLIPSE*/
		x = xU-xg[k];
		y = yU-yg[k];
                double EU = pow(b,2)*(pow(y*cos(theta[k]),2)+pow(x*sin(theta[k]),2)-y*x*sin(2*theta[k]))
                + pow(a,2)*(pow(x*cos(theta[k]),2)+pow(y*sin(theta[k]),2)+y*x*sin(2*theta[k]))
                - pow(a,2)*pow(b,2);

		x = xV-xg[k];
		y = yV-yg[k];
                double EV = pow(b,2)*(pow(y*cos(theta[k]),2)+pow(x*sin(theta[k]),2)-y*x*sin(2*theta[k]))
                + pow(a,2)*(pow(x*cos(theta[k]),2)+pow(y*sin(theta[k]),2)+y*x*sin(2*theta[k]))
                - pow(a,2)*pow(b,2);

		x = xS-xg[k];
		y = yS-yg[k];
                double ES = pow(b,2)*(pow(y*cos(theta[k]),2)+pow(x*sin(theta[k]),2)-y*x*sin(2*theta[k]))
                + pow(a,2)*(pow(x*cos(theta[k]),2)+pow(y*sin(theta[k]),2)+y*x*sin(2*theta[k]))
                - pow(a,2)*pow(b,2);

                Ip_S[k][i][j]=(ES <= 0);
                Ip_U[k][i][j]=(EU <= 0);
                Ip_V[k][i][j]=(EV <= 0);

		coloring[i][j] +=Ip_S[k][i][j];
#endif
#ifdef DISK
                chi_S[k][i][j]=((xS-xg[k])*(xS-xg[k])+(yS-yg[k])*(yS-yg[k])<= (rp[k]+data->eps)*(rp[k]+data->eps));
                chi_U[k][i][j]=((xU-xg[k])*(xU-xg[k])+(yU-yg[k])*(yU-yg[k])<= (rp[k]+data->eps)*(rp[k]+data->eps));
                chi_V[k][i][j]=((xV-xg[k])*(xV-xg[k])+(yV-yg[k])*(yV-yg[k])<= (rp[k]+data->eps)*(rp[k]+data->eps));

#ifndef SMOOTHING
                Ip_S[k][i][j]=((xS-xg[k])*(xS-xg[k])+(yS-yg[k])*(yS-yg[k])<= rp[k]*rp[k]);
                Ip_U[k][i][j]=((xU-xg[k])*(xU-xg[k])+(yU-yg[k])*(yU-yg[k])<= rp[k]*rp[k]);
                Ip_V[k][i][j]=((xV-xg[k])*(xV-xg[k])+(yV-yg[k])*(yV-yg[k])<= rp[k]*rp[k]);
#endif

#ifdef SMOOTHING

                //Smoothing S
                d = rp[k] - sqrt((xS-xg[k])*(xS-xg[k])+(yS-yg[k])*(yS-yg[k]));
                if( d < - data->eps)
                    Ip_S[k][i][j] = 0;
                else if( fabs(d) <= data->eps)
                    Ip_S[k][i][j] = .5*(1 + d/data->eps + (1./M_PI)*sin( M_PI* d/data->eps) );
                else if( d > data->eps)
                    Ip_S[k][i][j] = 1;

                //Smoothing U
                d = rp[k] - sqrt((xU-xg[k])*(xU-xg[k])+(yU-yg[k])*(yU-yg[k]));
                if( d < - data->eps)
                    Ip_U[k][i][j] = 0;
                else if( fabs(d) <= data->eps)
                    Ip_U[k][i][j] = .5*(1 + d/data->eps + (1./M_PI)*sin( M_PI* d/data->eps) );
                else if( d > data->eps)
                    Ip_U[k][i][j] = 1;

                //Smoothing V
                d = rp[k] - sqrt((xV-xg[k])*(xV-xg[k])+(yV-yg[k])*(yV-yg[k]));
                if( d < - data->eps)
                    Ip_V[k][i][j] = 0;
                else if( fabs(d) <= data->eps)
                    Ip_V[k][i][j] = .5*(1 + d/data->eps + (1./M_PI)*sin( M_PI* d/data->eps) );
                else if( d > data->eps)
                    Ip_V[k][i][j] = 1;

#endif

                xloc = xS-xg[k];
                yloc = yS-yg[k];
                delta = atan2(yloc, xloc); 
                coloring[i][j] += Ip_S[k][i][j];

                if((int) floor((delta-theta[k])/(M_PI/2.)) % 2 == 0 ){
                    coloring[i][j] = -coloring[i][j];
                }
#endif
                I_S[i][j] += Ip_S[k][i][j];
                I_U[i][j] += Ip_U[k][i][j];
                I_V[i][j] += Ip_V[k][i][j];
            }

            if(I_S[i][j] > 1 || I_U[i][j] > 1 || I_V[i][j] > 1 ){
                PetscPrintf(PETSC_COMM_WORLD, "Erreur : collision de particules \n");
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
    double* xg = data->xg;
    double* yg = data->yg;

    double*** chi_U = data->chi_U;
    double*** chi_V = data->chi_V;

    int m = data->m;
    int n = data->n;
    int Np = data->Np;
    double h = data->h;

    double xV, yU;
    for(int i=0; i<m; i++){
        xV = (i-0.5)*h;
        for(int j=0; j<n; j++){
            yU = (j-0.5)*h;
            u_s[i][j] = 0.;
            v_s[i][j] = 0.;
            for (int k = 0; k<Np; k++){
                u_s[i][j]+= chi_U[k][i][j]*(Up[k][3] - Omega_p[k][3]*(yU-yg[k]));
                v_s[i][j]+= chi_V[k][i][j]*(Vp[k][3] + Omega_p[k][3]*(xV-xg[k]));
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
    double** u_star = data->u_star;
    double** u_s = data->u_s;

    double** v_n = data->v_n;
    double** v_star = data->v_star;
    double** v_s = data->v_s;

    double** P = data->P;

    double Um = data->u_m;

    double H_U, H_U_old;
    double H_V, H_V_old;
    double lapU, lapV;
    double dpdx, dpdy;

    double uL, uR, vT, vB;
    double dudxR, dudxL, dudyT, dudyB, dvdxR, dvdxL, dvdyT, dvdyB;

    /* u_star ADAMS-BASHFORTH 2 */
    for (i=1; i<m-2; i++){
        for (j=1; j<n-1; j++){

            uR = .5*(u_n[i][j] + u_n[i+1][j]);
            uL = .5*(u_n[i][j] + u_n[i-1][j]);
            dudxR = (u_n[i+1][j]- u_n[i][j])/h;
            dudxL = (u_n[i][j]- u_n[i-1][j])/h;

            vT = .5*(v_n[i][j] + v_n[i+1][j]);
            vB = .5*(v_n[i][j-1] + v_n[i+1][j-1]);
            dudyT = (u_n[i][j+1] - u_n[i][j])/h;
            dudyB = (u_n[i][j] - u_n[i][j-1])/h;

            H_U = .5*(uR*dudxR + uL*dudxL) + .5*(vT*dudyT + vB*dudyB);

            if (data->iter == 1){
                H_U_old = H_U;
            }
            else
            {
                H_U_old = data->H_u_n_1[i][j];
            }

            // LAPLACIAN
            lapU = (u_n[i+1][j]+u_n[i-1][j]+u_n[i][j+1]+u_n[i][j-1]-4.*u_n[i][j])/(h*h);

            // PRESSURE term
            dpdx = (P[i+1][j]-P[i][j])/h;

#ifdef EXPLICIT
            // EXPLICIT VERSION
            u_star[i][j] = u_n[i][j] + dt*(-1.5*H_U + 0.5*H_U_old - dpdx + nu*lapU) - ramp*I_U[i][j]*(dt/dtau)*(u_n[i][j] - u_s[i][j]);
#else
            // IMPLICIT VERSION
            u_star[i][j] = (u_n[i][j] + dt*(-1.5*H_U + 0.5*H_U_old - dpdx + nu*lapU) + (dt/dtau)*ramp*I_U[i][j]*u_s[i][j])/(1.+ramp*I_U[i][j]*dt/dtau);
#endif
            data->H_u_n_1[i][j] = H_U;
        }
    }

    /*Outflow condition */
    for (j=1; j<n-1; j++){
        u_star[m-2][j] = u_n[m-2][j] - dt*Um*(u_n[m-2][j]-u_n[m-3][j])/h;
    }
    /*u_star[0][j] (inflow) is fixed once for all at the beginning */

    /* v_star  ADAMS-BASHFORTH 2 */
    for (i=1; i<m-1; i++){
        for (j=1; j<n-2; j++){

            uR = .5*(u_n[i][j] + u_n[i][j+1]);
            uL = .5*(u_n[i-1][j] + u_n[i-1][j+1]);
            dvdxR = (v_n[i+1][j]- v_n[i][j])/h;
            dvdxL = (v_n[i][j]- v_n[i-1][j])/h;

            vT = .5*(v_n[i][j] + v_n[i][j+1]);
            vB = .5*(v_n[i][j] + v_n[i][j-1]);
            dvdyT = (v_n[i][j+1] - v_n[i][j])/h;
            dvdyB = (v_n[i][j] - v_n[i][j-1])/h;

            H_V = .5*(uR*dvdxR + uL*dvdxL) + .5*(vT*dvdyT + vB*dvdyB);

            if (data->iter == 1){
                H_V_old = H_V;
            }
            else
            {
                H_V_old = data->H_v_n_1[i][j];
            }

            // LAPLACIAN
            lapV = (v_n[i+1][j]+v_n[i-1][j]+v_n[i][j+1]+v_n[i][j-1]-4.*v_n[i][j])/(h*h);

            // PRESSURE TERM
            dpdy = (P[i][j+1]-P[i][j])/h;

#ifdef EXPLICIT
            // EXPLICIT VERSION
            v_star[i][j] = v_n[i][j] + dt*(-1.5*H_V + 0.5*H_V_old - dpdy + nu*lapV) - ramp*I_V[i][j]*(dt/dtau)*(v_n[i][j] - v_s[i][j]);
#else
            // IMPLICIT VERSION
            v_star[i][j] = (v_n[i][j] + dt*(-1.5*H_V + 0.5*H_V_old - dpdy + nu*lapV) + (dt/dtau)*ramp*I_V[i][j]*v_s[i][j])/(1.+ramp*I_V[i][j]*dt/dtau);
#endif
            /* the value of v_star on the boundaries (j=0, j=n-2) is set to zero at allocation */

            data->H_v_n_1[i][j] = H_V;
        }
    }
}

void update_flow(Data* data) {

    int m = data->m;
    int n = data->n;
    double dt = data->dt;
    double h = data->h;

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
    //double y_ch, u_poiseuille;
    /* Set boundary for u_new, v_new */

    for (j = 1; j < n - 1; j++) {
        u_new[0][j] = u_n[0][j];
    }
    for (i = 1; i < m - 1; i++) {
        for (j = 1; j < n - 1; j++) {
            u_new[i][j] = u_star[i][j] - dt * (phi[i + 1][j] - phi[i][j]) / h;
            v_new[i][j] = v_star[i][j] - dt * (phi[i][j + 1] - phi[i][j]) / h;
            P[i][j] += phi[i][j];
        }
    }

#ifdef TEMP
    update_scalars(data);
#endif

    free2Darray(u_n_1, m);
    free2Darray(v_n_1, m);

    data->u_n_1 = u_n;
    data->u_n = u_new;

    data->v_n_1 = v_n;
    data->v_n = v_new;

}

void update_scalars(Data* data)
{
    /* TEMPERATURE AND SPECIES */

    int m = data->m;
    int n = data->n;
    double h = data->h;
    int i,j;

    double **T_new = make2DDoubleArray(m, n);
    double **T_n = data->T_n;
    double **T_n_1 = data->T_n_1;
    double **Ts = data->Ts;

    double **u_n = data->u_n;
    double **v_n = data->v_n;

    int Ns = data->Ns;
    double ***C_new = make3DDoubleArray(Ns, m, n);
    double ***C_n = data->C_n;
    double ***C_n_1 = data->C_n_1;
    double ***Cs = data->Cs;

    double alpha_f = data->alpha_f;
    double *Df = data->Df;
    double **I_S = data->I_S;
    double ramp = data->ramp;
    double dtau = data->dtau;
    double dt = data->dt;

    double H_T, H_T_old, lapT;
    double H_C, H_C_old, lapC;

    /** TEMPERATURE **/

    for (i = 1; i < m - 1; i++) {
        for (j = 1; j < n - 1; j++) {

            // ADVECTIVE TERMS
            H_T = .5*(u_n[i][j]*(T_n[i+1][j]-T_n[i][j])/h + u_n[i-1][j]*(T_n[i][j]-T_n[i-1][j])/h)
                  + .5*(v_n[i][j]*(T_n[i][j+1]-T_n[i][j])/h + v_n[i][j-1]*(T_n[i][j]-T_n[i][j-1])/h);

            if(data->iter == 1)
            {
                H_T_old = H_T;
            }
            else{
                H_T_old = data->H_T_n_1[i][j];
            }

            // DIFFUSION TERM
            lapT = (T_n[i + 1][j] + T_n[i - 1][j] + T_n[i][j + 1] + T_n[i][j - 1] - 4. * T_n[i][j]) / (h * h);

#ifdef EXPLICIT
            // EXPLICIT VERSION
            T_new[i][j] = T_n[i][j] + dt * (-1.5 * H_T + 0.5 * H_T_old
                                            + alpha_f * lapT)
                                    - ramp*I_S[i][j]*(dt/dtau)*(T_n[i][j] - Ts[i][j]);
#else
            // IMPLICIT VERSION
            T_new[i][j] = (T_n[i][j] + dt * (-1.5 * H_T + 0.5 * H_T_old
                                             + alpha_f * lapT
                                             + ramp * I_S[i][j] * Ts[i][j] / dtau)) /
                          (1. + ramp * I_S[i][j] * dt / dtau);
#endif
            data->H_T_n_1[i][j] = H_T;

            /** SPECIES **/

            for (int s = 0; s < Ns-1; s++) {
                // ADVECTIVE TERMS
                H_C = .5*(u_n[i][j]*(C_n[s][i+1][j]-C_n[s][i][j])/h + u_n[i-1][j]*(C_n[s][i][j]-C_n[s][i-1][j])/h)
                      + .5*(v_n[i][j]*(C_n[s][i][j+1]-C_n[s][i][j])/h + v_n[i][j-1]*(C_n[s][i][j]-C_n[s][i][j-1])/h);

                if(data->iter == 1)
                {
                    H_C_old = H_C;
                }
                else{
                    H_C_old = data->H_C_n_1[s][i][j];
                }

                // DIFFUSION TERM
                lapC = (C_n[s][i + 1][j] + C_n[s][i - 1][j] + C_n[s][i][j + 1] + C_n[s][i][j - 1] - 4. * C_n[s][i][j]) / (h * h);
#ifdef EXPLICIT
              //EXPLICIT VERSION
              C_new[s][i][j] = C_n[s][i][j] + dt * (-1.5 * H_C + 0.5 * H_C_old
                                                     + Df[s] * lapC)
                                              - ramp*I_S[i][j]*(dt/dtau)*(C_n[s][i][j] - Cs[s][i][j]);
#else
                // IMPLICIT VERSION
                C_new[s][i][j] = (C_n[s][i][j] + dt * (-1.5 * H_C + 0.5 * H_C_old
                                                       + Df[s] * lapC
                                                       + ramp * I_S[i][j] * Cs[s][i][j] / dtau))/
                                 (1. + ramp * I_S[i][j] * dt / dtau);

#endif
                data->H_C_n_1[s][i][j] = H_C;
            }

        }
    }

    free2Darray(T_n_1, m);
    free3Darray(C_n_1, Ns, m);

    data->T_n_1 = T_n;
    data->T_n = T_new;

    data->C_n_1 = C_n;
    data->C_n = C_new;
}

void get_vorticity(Data* data){
    int m = data->m;
    int n = data->n;
    double h = data->h;

    int i,j;

    double** omega = data->omega;
    double** u_n = data->u_n;
    double** v_n = data->v_n;

    for (i = 0; i < m - 1; i++) {
        for (j = 1; j < n - 2; j++) {
            omega[i][j] = (v_n[i + 1][j] - v_n[i][j]) / h - (u_n[i][j + 1] - u_n[i][j]) / h;
        }
    }
}

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
