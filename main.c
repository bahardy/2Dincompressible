#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <petsc.h>
#include <sys/stat.h>
#include "main.h"
#include "write.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

//#define RECOVER
//#define MOVE
#define TEMP
#define TWO_WAY
#define RAMPING
#define WRITE
#define DISK
#define SLIP
#define CHANNEL
//#define GRAVITY
//#define SMOOTHING
//#define ELLIPSE


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

    /* DIMENSIONS */
    data.Dp = 1.;
    data.d = 5.*data.Dp;
    data.H = 0.5*data.d;
    data.L = 5.*data.Dp;
    data.h = data.Dp/40;
    data.eps = 0;
#ifdef SMOOTHING
    data.eps = data.h;
#endif
    /* NON-DIMENSIONAL NUMBERS */
    data.Pr = 0.7;
    data.Le = 1; /* Lewis number, ratio between Sc and Prandtl */
    data.Sc = data.Le*data.Pr;
    data.Rep = 40.;
    data.Fr = sqrt(1e3);

    /* FLOW */
    data.u_m = 1.;
    //data.u_max = 1.5*data.u_m;
    //data.nu = 0.03926;
    data.nu = data.u_m*data.Dp/data.Rep;
    data.g = 0;
#ifdef GRAVITY
    data.g = 9.81;
#endif
    /* ENERGY */
    data.alpha_f = data.nu/data.Pr;
    data.Tm0 = 1; // cup-mixing temperature at the inlet
    data.Tp0 = 0.5; // initial particle temperature

    /* PHYSICAL PARAMETERS */
    data.rho_f = 1.;
    data.rho_p = 1000.;
    data.rho_r = data.rho_p/data.rho_f;
    data.cp = 1000.;
    data.cf = 1000.;
    data.cr = data.cp/data.cf;

    /* SPECIES */
    data.Ns = 2;
    data.Np = 1;
    data.Df = make1DDoubleArray(data.Ns);
    data.Df[0] = data.nu/data.Sc;
    data.Df[1] = data.nu/data.Sc;
    data.dH = 0;
    data.CA0 = 1.;
    data.CB0 = 0.;

    /* GRID */
    data.N = (int) (data.d/data.h);
    data.M = (int) (data.L/data.h);
    data.n = data.N + 2; /*for ghost points */
    data.m = data.M + 2; /*for ghost points */


    /* TIME INTEGRATION */
    data.CFL = 0.5; /*Courant-Freidrichs-Lewy condition on convective term */
    data.r = .25; /* Fourier condition on diffusive term */
    double dt_CFL = data.CFL*data.h/data.u_m;
    double dt_diff = data.r*data.h*data.h/data.nu;

    data.ratio_dtau_dt = 1e-3;
    data.dt = fmin(dt_CFL, dt_diff);
    data.dtau = data.ratio_dtau_dt*data.dt;

    if(rank == 0){
        printf("Rep = %f\n", data.Rep);
        printf("ratio L/d = %f \n", data.L/data.d);
        printf("Um = %f\n", data.u_m);
        printf("dt_CFL = %f\n", dt_CFL);
        printf("dt_diff = %f\n", dt_diff);
        printf("CFL condition set to %f h/Umax \n", data.CFL);
        printf("dt = %f\n", data.dt);
        printf("dtau = %f\n", data.dtau);
        printf("ratio dtau/dt = %f \n", data.ratio_dtau_dt);
    }

    /* Writing */

    if (argc >= 3) {
        sscanf(argv[1], "%d", &(data.T_write));
        sscanf(argv[2], "%d", &(data.N_write));
    }
    else{
        data.T_write = 1; /* number of time steps between two writings */
        data.N_write = 50; /* number of times we write in files */
    }

    double Tf = data.N_write*data.T_write*data.dt;
    data.Tf = Tf;
    data.t_move = 0; //data.Tf/10.;
    data.nKmax = 2;
    data.Kmax = 50; /* number of ramping steps */

    double t, t_start;
    int iter_start;
    data.ramp = 1;
    double surf = 0.;

    if(rank == 0){
        printf("Write every %d * dt \n", data.T_write);
        printf("Write %d times \n", data.N_write);
        printf("Final time : %f \n \n", data.Tf);
    }

    /**------------------------------- Matrix Fields Creation  ------------------------------- **/

    int m = data.m;
    int n = data.n;
    int Np = data.Np;
    int Ns = data.Ns;

    data.coloring = make2DDoubleArray(m,n);
    data.C_n = make3DDoubleArray(Ns,m,n);
    data.C0 = make1DDoubleArray(Ns); // inlet concentration
    data.C_n_1 = make3DDoubleArray(Ns,m,n);
    data.Cs = make3DDoubleArray(Ns,m,n);
    data.Cp = make2DDoubleArray(Np,Ns);
    data.dudt = make1DDoubleArray(Np);
    data.dvdt = make1DDoubleArray(Np);
    data.domegadt = make1DDoubleArray(Np);
    data.dTdt = make1DDoubleArray(Np);
    data.dCdt = make2DDoubleArray(Np,Ns);
    data.dp = make1DDoubleArray(Np);
    data.F = make2DDoubleArray(Np,3);
    data.Fx = make1DDoubleArray(Np);
    data.Fy = make1DDoubleArray(Np);
    data.G = make2DDoubleArray(Np,3);
    data.Ip_S = make3DDoubleArray(Np,m,n);
    data.Ip_U = make3DDoubleArray(Np,m,n);
    data.Ip_V = make3DDoubleArray(Np,m,n);
    data.I_S = make2DDoubleArray(m,n);
    data.I_U = make2DDoubleArray(m,n);
    data.I_V = make2DDoubleArray(m,n);
    data.II = make1DDoubleArray(Np);
    data.J = make1DDoubleArray(Np);
    data.Mz = make2DDoubleArray(Np,3);
    data.omega = make2DDoubleArray(m,n);
    data.Omega_p = make2DDoubleArray(Np,4);
    data.phi = make2DDoubleArray(m,n);
    data.Qm = make2DDoubleArray(Np,Ns);
    data.P = make2DDoubleArray(m,n);
    data.PP = make3DDoubleArray(Np, Ns, 3);
    data.Q = make1DDoubleArray(Np);
    data.QQ = make2DDoubleArray(Np,3);
    data.Qr = make2DDoubleArray(Np,3);
    data.rp = make1DDoubleArray(Np);
    data.Reh = make2DDoubleArray(m,n);
    data.Reh_omega = make2DDoubleArray(m,n);
    data.Sp = make1DDoubleArray(Np);
    data.theta = make1DDoubleArray(Np);

    data.T_n = make2DDoubleArray(m,n);
    data.T_n_1 = make2DDoubleArray(m,n);
    data.Tp = make1DDoubleArray(Np);
    data.Ts = make2DDoubleArray(m,n);
    data.Tz = make1DDoubleArray(Np);

    data.u_n = make2DDoubleArray(m,n);
    data.u_n_1 = make2DDoubleArray(m,n);
    data.Up = make2DDoubleArray(Np,4);
    data.u_s = make2DDoubleArray(m,n);
    data.u_star = make2DDoubleArray(m,n);

    data.v_n = make2DDoubleArray(m,n);
    data.v_n_1 = make2DDoubleArray(m,n);
    data.Vp = make2DDoubleArray(Np,4);
    data.v_s = make2DDoubleArray(m,n);
    data.v_star = make2DDoubleArray(m,n);

    data.xg = make1DDoubleArray(Np);
    data.yg = make1DDoubleArray(Np);


    /** ------------------------------- Fields Initialization ------------------------------- **/

    /* Particles position */
    data.xg[0] = data.d;
    data.yg[0] = data.H;
    data.dp[0] = data.Dp;
    data.rp[0] = .5*data.Dp;
    data.theta[0] = 0; // M_PI/10.

#ifndef TWO_WAY
    //impulsively started cylinder : we impose the motion
    data.Up[0][0] = data.u_m;
    data.Up[0][1] = data.Up[0][0];
    data.Up[0][2] = data.Up[0][0];
    data.Up[0][3] = data.Up[0][2];
#endif

    for(int k=0; k<Np; k++){
#ifdef DISK
        data.Sp[k]=M_PI*data.rp[k]*data.rp[k];
        data.II[k]=(data.dp[k]*data.dp[k])/8.; /* IN 2-D !! */
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

    FILE* fichier_stat = fopen("results/stats.txt", "w+");
    FILE* fichier_forces = fopen("results/forces.txt", "w+");
    FILE* fichier_fluxes = fopen("results/fluxes.txt", "w+");
    FILE* fichier_particle = fopen("results/particle.txt", "w+");

    /** ------------------------------- INITIALIZATION of the domain ------------------------------- **/

    //double y_ch;
    /* VELOCITY : horizontal flow Um  */
    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
//            y_ch = (j-0.5)*data.h - data.H;
//            data.u_n[i][j] = data.u_max*(1.-(y_ch/data.H)*(y_ch/data.H));
            data.u_n[i][j] = data.u_m;
            data.u_n_1[i][j] = data.u_n[i][j];
            data.u_star[i][j] = data.u_n[i][j];
            data.T_n[i][j] = data.Tm0;
            data.T_n_1[i][j] = data.T_n[i][j];
            /* v_n is initially at zero */
        }
    }

    /** ----- BOUNDARY CONDITION -----------**/

    //y_ch = (j-0.5)*data.h - data.H;
    //data.u_max*(1.-(y_ch/data.H)*(y_ch/data.H));

    // Ghost point for u_n_1 (first time EE)
    for (int i=0; i<m; i++) {
#ifdef SLIP
	// slip walls : dudn = 0
        data.u_n_1[i][0] = data.u_n_1[i][1];
        data.u_n_1[i][n-1] = data.u_n_1[i][n-2];
#endif
#ifndef SLIP
	// no-slip walls : u = 0
        data.u_n_1[i][0] = -0.2 *(data.u_n_1[i][3] - 5. * data.u_n_1[i][2] + 15. * data.u_n_1[i][1]);
        data.u_n_1[i][n - 1] = -0.2 * (data.u_n_1[i][n - 4] - 5. * data.u_n_1[i][n - 3] + 15. * data.u_n_1[i][n - 2]);
#endif
        // Ghost point for T_n_1 (first time EE)
        data.T_n_1[i][0] = data.T_n_1[i][1];
        data.T_n_1[i][n-1] = data.T_n_1[i][n-2];
    }

    /*Initialization of particles temperatures */
    for(int k=0; k<Np; k++){
        data.Tp[k] = data.Tp0;
    }

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

        for (int k = 0; k < Np; k++) {
            get_masks(&data);
            integrate_penalization(&data, &surf, k);
            /* dudt, dvdt, etc = 0 because the particle is fixed */
            compute_forces_fluxes(&data, k);
        }

        get_ghosts(&data, data.Tm0, data.C0);
        get_Ustar_Vstar(&data, data.ramp);
        poisson_solver(&data, rank, nbproc);
        update_flow(&data);
        diagnostic(&data);

    }
#endif

#ifdef WRITE
    /*INITIAL SOLUTION (t=0) AFTER RAMPING */
    if(rank==0){
        writeFields(&data, 0);
	    //writeMask(&data);
    }
#endif


    /** -------------------------------TIME STEPPING FROM BEGINNING ------------------------------- **/
    iter_start = 1;
    t_start = 0.;

    /** -------------------------------TIME STEPPING ------------------------------- **/
    data.iter = iter_start;
    t = t_start;

    // we feed reactants at the inlet
    data.C0[0] = data.CA0;
    data.C0[1] = data.CB0;

    while(t < data.Tf){

        PetscPrintf(PETSC_COMM_WORLD, "\n \n BEGIN iter %d : t = %f \n", data.iter, t);

        /** --- SOLVE PARTICULAR PHASE --- */
        int flag_out = 0;
        int k;
        for (k = 0; k<Np; k++){
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
            update_Tp(&data, k);
            update_Cp(&data, k);
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
        get_Cs(&data);
#endif
        get_Ustar_Vstar(&data, data.ramp);
        clock_t t_init = clock();
        poisson_solver(&data, rank, nbproc);
        clock_t t_final = clock();
        double t_Poisson = ((double) (t_final - t_init))/CLOCKS_PER_SEC;
        PetscPrintf(PETSC_COMM_WORLD, "Poisson solver took %f seconds \n", t_Poisson);

        update_flow(&data);
        get_ghosts(&data, data.Tm0, data.C0);
        get_vorticity(&data);

        diagnostic(&data);


#ifdef WRITE
        if(rank == 0){
            fprintf(fichier_stat, "%3.13e \t  %3.13e \n", data.Reh_max, data.Reh_omega_max);
            fflush(fichier_stat);
            fprintf(fichier_forces, "%3.13e \t  %3.13e  \t %3.13e \n", data.Fx[0],  data.Fy[0], data.Tz[0]);
            fflush(fichier_forces);
            fprintf(fichier_fluxes, "%3.13e \t  %3.13e \n", data.Q[0], data.Qm[0][0]);
            fflush(fichier_fluxes);
            fprintf(fichier_particle, "%3.13e \t  %3.13e \t %3.13e \t %3.13e \t %3.13e  \t %3.13e \t %3.13e \n", data.xg[0], data.yg[0], data.theta[0],
                    data.Up[0][3], data.Vp[0][3], data.Omega_p[0][3], data.Tp[0]);
            fflush(fichier_particle);

            if(data.iter % data.T_write == 0){
                writeFields(&data, data.iter);
            }
        }
#endif
        /* Increment time */
        t += data.dt;
        data.iter ++;

    }

    fclose(fichier_forces);
    fclose(fichier_fluxes);
    fclose(fichier_particle);
    fclose(fichier_stat);

    /* Free memory */
    free(data.xg), free(data.yg), free(data.theta), free(data.dp), free(data.rp), free(data.Sp), free(data.II), free(data.J);
    free(data.dudt), free(data.dvdt), free(data.domegadt), free(data.dTdt); free2Darray(data.dCdt, Np);
    free(data.Fx), free(data.Fy), free(data.Tz), free(data.Q), free2Darray(data.Qm, Np);
    free2Darray(data.u_n,m), free2Darray(data.u_n_1,m), free2Darray(data.u_star,m), free2Darray(data.u_s,m);
    free2Darray(data.v_n,m), free2Darray(data.v_n_1,m), free2Darray(data.v_star,m), free2Darray(data.v_s,m);
    free2Darray(data.omega, m); free2Darray(data.Reh,m); free2Darray(data.Reh_omega,m);
    free2Darray(data.P,m), free2Darray(data.phi, m);
    free2Darray(data.T_n,m),  free2Darray(data.T_n_1,m), free2Darray(data.Ts, m);
    free3Darray(data.C_n, Ns, m), free3Darray(data.C_n_1, Ns, m), free3Darray(data.Cs, Ns, m), free(data.C0);
    free2Darray(data.Up,Np), free2Darray(data.Vp, Np), free2Darray(data.Omega_p,Np), free(data.Tp), free2Darray(data.Cp, Np);
    free2Darray(data.F, Np), free2Darray(data.G, Np), free2Darray(data.Mz, Np), free2Darray(data.QQ, Np), free3Darray(data.PP, Np,Ns), free2Darray(data.Qr, Np);
    free2Darray(data.I_S, m), free2Darray(data.I_U, m), free2Darray(data.I_V, m), free2Darray(data.coloring, m);
    free3Darray(data.Ip_S, Np,m), free3Darray(data.Ip_U, Np,m), free3Darray(data.Ip_V, Np, m);

    PetscFinalize();
    return 0;
}

void compute_forces_fluxes(Data* data, int k)
{

    double* Fx = data->Fx;
    double* Fy = data->Fy;
    double* Tz = data->Tz;
    double* Q = data->Q;
    double** Qm = data->Qm;

    double** F = data->F;
    double** G = data->G;
    double** M = data->Mz;
    double** QQ = data->QQ;
    double*** PP = data->PP;

    double* dudt = data->dudt;
    double* dvdt = data->dvdt;
    double* domegadt = data->domegadt;
    double* dTdt = data->dTdt;

    double* Sp = data->Sp;
    double* II = data->II;

    double rho_f = data->rho_f;
    double cf = data->cf;

    PetscPrintf(PETSC_COMM_WORLD,"\n F integration = %1.6e \n", k+1, F[k][2]);
    PetscPrintf(PETSC_COMM_WORLD,"G integration = %1.6e \n", k+1, G[k][2]);
    PetscPrintf(PETSC_COMM_WORLD,"M integration = %1.6e \n", k+1, M[k][2]);

    PetscPrintf(PETSC_COMM_WORLD,"dudt = %1.6e \n", k+1, dudt[k]);
    PetscPrintf(PETSC_COMM_WORLD,"dvdt = %1.6e \n", k+1, dvdt[k]);
    PetscPrintf(PETSC_COMM_WORLD,"domegadt = %1.6e \n", k+1, domegadt[k]);

    Fx[k] = rho_f*(Sp[k]*dudt[k] + F[k][2]);
    Fy[k] = rho_f*(Sp[k]*dvdt[k] + G[k][2]);
    Tz[k] = rho_f*(Sp[k]*II[k]*domegadt[k] + M[k][2]);
    Q[k] = rho_f*cf*(Sp[k]*dTdt[k] + QQ[k][2]);
    Qm[k][0] = PP[k][0][2];


    PetscPrintf(PETSC_COMM_WORLD,"Hydrodynamic force along -x dir on particle %d = %1.6e [N/m]  \n", k+1, Fx[k]);
    PetscPrintf(PETSC_COMM_WORLD,"Hydrodynamic force along -y dir on particle %d = %1.6e [N/m]  \n", k+1, Fy[k]);
    PetscPrintf(PETSC_COMM_WORLD,"Torque on particle %d = %1.6e [N]  \n", k+1, Tz[k]);
    //PetscPrintf(PETSC_COMM_WORLD,"Heat flux on particle %d = %1.6e [W/m] \n", k+1, Q[k]);
    //PetscPrintf(PETSC_COMM_WORLD,"Molar flux of A on particle %d = %1.6e [mol/(m.s)] \n", k+1, Phi[k][0]);
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
    double h = data->h;
    double nu = data->nu;
    double m = data->m;
    double n = data->n;

    data->Reh_max = 0.;
    data->Reh_omega_max = 0.;

    int i, j;
    for (i=1; i<m-1; i++) {
        for (j = 1; j<n-1; j++) {
            Reh[i][j] = (fabs(u_n[i][j]) + fabs(v_n[i][j]))*h/nu;
            Reh_omega[i][j] = fabs(omega[i][j])*h*h/n;
            if(Reh[i][j] > data->Reh_max){
                data->Reh_max = Reh[i][j];
            }
            if(Reh_omega[i][j] > data->Reh_omega_max){
                data->Reh_omega_max = Reh_omega[i][j];
            }
        }
    }


}

int integrate_penalization(Data *data, double* surf, int k)
{
    // Integral terms
    double** F = data->F;
    double** G = data->G;
    double** M = data->Mz;
    double** QQ = data-> QQ;
    double*** PP = data-> PP;

    double*** Ip_U = data->Ip_U;
    double*** Ip_V = data->Ip_V;
    double*** Ip_S = data->Ip_S;
    double** u_n = data->u_n;
    double** v_n = data->v_n;
    double*** C_n = data->C_n;
    double** T_n = data->T_n;
    double** u_s = data->u_s;
    double** v_s = data->v_s;
    double** Ts = data->Ts;
    double*** Cs = data->Cs;
    double* xg = data->xg;
    double* yg = data->yg;
    double* rp = data->rp;
    int Ns = data->Ns;
    int m = data->m;
    int n = data->n;
    double h = data->h;
    double dtau = data->dtau;

    /* Force along x-direction */
    F[k][0] = F[k][1]; /* n-2*/
    F[k][1] = F[k][2]; /* n-1*/
    //F[k][2] = 0.; /* n*/

    /* Force along y-direction */
    G[k][0] = G[k][1];
    G[k][1] = G[k][2];
    //G[k][2] = 0.;

    /* Moment along z-direction */\
    M[k][0] = M[k][1];
    M[k][1] = M[k][2];
    //M[k][2] = 0.;

#ifdef TEMP
    /* Particle heat balance  */
    QQ[k][0] = QQ[k][1]; /* n-2*/
    QQ[k][1] = QQ[k][2]; /* n-1*/
    //QQ[k][2] = 0.; /* n*/

    for(int s=0; s<Ns; s++){
        PP[k][s][0] = PP[k][s][1];
        PP[k][s][1] = PP[k][s][2];
     //  PP[k][s][2] = 0.;
    }
#endif

    int startX = (int) floor((xg[k] - rp[k]) / h);
    PetscPrintf(PETSC_COMM_WORLD, "startX = %d \t", startX);
    int endX = (int) ceil((xg[k]+rp[k])/h);
    PetscPrintf(PETSC_COMM_WORLD,"endX = %d \t", endX);

    int startY = (int) floor((yg[k]-rp[k])/h);
    PetscPrintf(PETSC_COMM_WORLD,"startY = %d \t", startY);
    int endY = (int) ceil((yg[k]+rp[k])/h);
    PetscPrintf(PETSC_COMM_WORLD,"endY = %d \t \n", endY);

    if(startY <= 1||endY >= n-1){
        PetscPrintf(PETSC_COMM_WORLD,"Wall collision! \n");
        return 1;
    }

    if(endX >= m-1){
        PetscPrintf(PETSC_COMM_WORLD,"Particle leaves the channel! \n");
        return 1;
    }

    double Fint, Gint, Mint, Qint, *Qmint, sint;
    Fint = 0.;
    Gint = 0.;
    Mint = 0.;
    Qint = 0.;
    Qmint = make1DDoubleArray(Ns);
    sint = 0.;

    double h2;
    h2 = h*h;
    double yU, xV, f, g;

    //for(int i=startX; i<=endX; i++){//
    for(int i = 0; i<m; i++){
        xV = (i-0.5)*h;
        //for(int j=startY; j<=endY; j++){
        for(int j=0; j<n; j++){
            yU = (j-0.5)*h;
            f = -Ip_U[k][i][j]*(u_n[i][j]-u_s[i][j]);
            g = -Ip_V[k][i][j]*(v_n[i][j]-v_s[i][j]);
#ifdef TEMP
            double q = -Ip_S[k][i][j]*(T_n[i][j]-Ts[i][j]);
            double* qm = make1DDoubleArray(Ns);
            for(int s=0; s<Ns; s++){
                qm[s] = -Ip_S[k][i][j]*(C_n[s][i][j]-Cs[s][i][j]);
            }
#endif
            sint += Ip_S[k][i][j]*h*h;
            Fint += f; /* units : m/s */
            Gint += g; /* units : m/s */
            Mint += ((xV-xg[k])*g-(yU-yg[k])*f);/* units: m^2/s */
#ifdef TEMP
            Qint += q; /*units : K */
            for(int s=0; s<Ns; s++){
                Qmint[s] += qm[s]; /*units : mol/m.s */
            }
#endif
        }
    }
    Fint *= h2/dtau; /* units : m^3/s; */
    Gint *= h2/dtau;
    Mint *= h2/dtau;
    Qint *= h2/dtau; /* units : K*m^2/s */

    F[k][2] = -Fint;
    G[k][2] = -Gint;
    M[k][2] = -Mint;
    QQ[k][2] = -Qint;
    *surf = sint;

#ifdef TEMP
    for(int s=0; s<Ns; s++){
        PP[k][s][2] = -Qmint[s]*h2/dtau;
    }
#endif
    PetscPrintf(PETSC_COMM_WORLD, "Particle surface is %f\n", *surf);
    return 0;
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
#ifdef CAVITY
        /* slip */
        v_n[0][j] = v_n[1][j];
#endif
#ifdef CHANNEL
        /* Inflow : horizontal flow --> v_n = 0 */
        v_n[0][j] = -0.2*(v_n[3][j] - 5.*v_n[2][j] + 15.*v_n[1][j]);
#endif
        /* Natural outflow : dV/dx =0 */
        v_n[m-1][j] = v_n[m-2][j];

#ifdef TEMP
        /* On T_n and C */
        /* Inflow : T_n uniform  */
        T_n[0][j] = -0.2*(T_n[3][j]-5.*T_n[2][j]+15.*T_n[1][j]-16.*T0);

        /* Inflow : CA = CA0; CB = CB0 */
        C[0][0][j] = -0.2*(C[0][3][j]-5.*C[0][2][j]+15.*C[0][1][j]-16.*C0[0]);
        C[1][0][j] = -0.2*(C[1][3][j]-5.*C[1][2][j]+15.*C[1][1][j]-16.*C0[1]);

        /*Outflow : We cancel axial dispersion d2T/dx2 = 0; d2C/dx2 = 0; */
        T_n[m-1][j] = (7.*T_n[m-2][j]-5.*T_n[m-3][j]+T_n[m-4][j])/3.;
        for(int s=0; s<Ns; s++){
            C[s][m-1][j] = (7.*C[s][m-2][j]-5.*C[s][m-3][j]+C[s][m-4][j])/3.;
        }
#endif

    }
}

void get_masks(Data* data)
{
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
    //double dist;
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


#ifndef SMOOTHING
                Ip_S[k][i][j]=((xS-xg[k])*(xS-xg[k])+(yS-yg[k])*(yS-yg[k])<= rp[k]*rp[k]);
                Ip_U[k][i][j]=((xU-xg[k])*(xU-xg[k])+(yU-yg[k])*(yU-yg[k])<= rp[k]*rp[k]);
                Ip_V[k][i][j]=((xV-xg[k])*(xV-xg[k])+(yV-yg[k])*(yV-yg[k])<= rp[k]*rp[k]);
#endif

#ifdef SMOOTHING

                //Smoothing S
                dist = rp[k] - sqrt((xS-xg[k])*(xS-xg[k])+(yS-yg[k])*(yS-yg[k]));
                if( dist < - data->eps)
                    Ip_S[k][i][j] = 0;
                else if( fabs(dist) <= data->eps)
                    Ip_S[k][i][j] = .5*(1 + dist/data->eps + (1./M_PI)*sin( M_PI* dist/data->eps) );
                else if( dist > data->eps)
                    Ip_S[k][i][j] = 1;

                //Smoothing U
                dist = rp[k] - sqrt((xU-xg[k])*(xU-xg[k])+(yU-yg[k])*(yU-yg[k]));
                if( dist < - data->eps)
                    Ip_U[k][i][j] = 0;
                else if( fabs(dist) <=data->eps)
                    Ip_U[k][i][j] = .5*(1 + dist/data->eps + (1./M_PI)*sin( M_PI* dist/data->eps) );
                else if( dist > data->eps)
                    Ip_U[k][i][j] = 1;

                //Smoothing V
                dist = rp[k] - sqrt((xV-xg[k])*(xV-xg[k])+(yV-yg[k])*(yV-yg[k]));
                if( dist < - data->eps)
                    Ip_V[k][i][j] = 0;
                else if( fabs(dist) <= data->eps)
                    Ip_V[k][i][j] = .5*(1 + dist/data->eps + (1./M_PI)*sin( M_PI* dist/data->eps) );
                else if( dist > data->eps)
                    Ip_V[k][i][j] = 1;

#endif

                xloc = xS-xg[k];
                yloc = yS-yg[k];
                delta = atan2(yloc, xloc); 
                coloring[i][j] += Ip_S[k][i][j];

                if((int)((delta-theta[k])/(M_PI/2.)) % 2 == 0 ){
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
    double*** Ip_S = data->Ip_S;
    int m = data->m;
    int n = data->n;
    int Np = data->Np;

    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            Cs[1][i][j] = 0.;
            for(int k=0; k<Np; k++){
                Cs[1][i][j] += Ip_S[k][i][j]*Cp[k][1];
            }
        }
    }

}
void get_Ts(Data* data)
{
    double** Ts = data-> Ts;
    double* Tp = data->Tp;
    double*** Ip_S = data->Ip_S;
    int m = data->m;
    int n = data->n;
    int Np = data->Np;

    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            Ts[i][j] = 0.;
            for(int k = 0; k<Np; k++){
                Ts[i][j]+= Ip_S[k][i][j]*Tp[k];

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

    double*** Ip_U = data->Ip_U;
    double*** Ip_V = data->Ip_V;

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
                u_s[i][j]+= Ip_U[k][i][j]*(Up[k][3] - Omega_p[k][3]*(yU-yg[k]));
                v_s[i][j]+= Ip_V[k][i][j]*(Vp[k][3] + Omega_p[k][3]*(xV-xg[k]));
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
    //double H = data->H;
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

    double Um = data->u_m;

    double Uij, Vij, Uij_old, Vij_old;
    double H_U, H_U_old;
    double H_V, H_V_old;
    double lapU, lapV;
    double dpdx, dpdy;
    //double Upoiseuille, y_ch;


    /* u_star ADAMS-BASHFORTH 2 */
    for (i=1; i<m-2; i++){
        for (j=1; j<n-1; j++){
            //only interior points
            // average vertical velocity
            Vij = .25*(v_n[i][j]+v_n[i+1][j]+v_n[i][j-1]+v_n[i+1][j-1]); /* average of the four neighbour points */
            Vij_old = .25*(v_n_1[i][j]+v_n_1[i+1][j]+v_n_1[i][j-1]+v_n_1[i+1][j-1]);
            // advective terms
            H_U =  u_n[i][j]*(u_n[i+1][j]-u_n[i-1][j])/(2.*h) + Vij*(u_n[i][j+1]-u_n[i][j-1])/(2.*h);
            H_U_old = u_n_1[i][j]*(u_n_1[i+1][j]-u_n_1[i-1][j])/(2.*h) + Vij_old*(u_n_1[i][j+1]-u_n_1[i][j-1])/(2.*h);
            // laplacian
            lapU = (u_n[i+1][j]+u_n[i-1][j]+u_n[i][j+1]+u_n[i][j-1]-4.*u_n[i][j])/(h*h);
            // pressure term
            dpdx = (P[i+1][j]-P[i][j])/h;

            u_star[i][j] = (u_n[i][j] + dt*(-1.5*H_U + 0.5*H_U_old - dpdx + nu*lapU) + (dt/dtau)*ramp*I_U[i][j]*u_s[i][j])/(1.+ramp*I_U[i][j]*dt/dtau);
        }
    }

#ifdef CAVITY
    for (j=1; j<n-1; j++){
        u_star[0][j] = u_n[1][j];//(4.*u_n[1][j] - u_n[2][j])/3.;
        u_star[m-2][j] = u_n[m-3][j];//(4.*u_n[m-3][j] - u_n[m-4][j])/3.;
    }
#endif

#ifdef CHANNEL
    /*Outflow condition */
    for (j=1; j<n-1; j++){
        u_star[m-2][j] = u_n[m-2][j] - dt*Um*(u_n[m-2][j]-u_n[m-3][j])/h;
    }
#endif
    /*u_star[0][j] (inflow) is fixed once for all at the beginning */

    /* v_star  ADAMS-BASHFORTH 2 */
    for (i=1; i<m-1; i++){
        for (j=1; j<n-2; j++){
            Uij = .25*(u_n[i][j]+u_n[i-1][j]+u_n[i][j+1]+u_n[i-1][j+1]); /* average of the four neighbour points */
            Uij_old = .25*(u_n_1[i][j]+u_n_1[i-1][j]+u_n_1[i][j+1]+u_n_1[i-1][j+1]);
            //Advective terms
            H_V = Uij*(v_n[i+1][j]-v_n[i-1][j])/(2.*h) + v_n[i][j]*(v_n[i][j+1]-v_n[i][j-1])/(2.*h);
            H_V_old = Uij_old*(v_n_1[i+1][j]-v_n_1[i-1][j])/(2.*h) + v_n_1[i][j]*(v_n_1[i][j+1]-v_n_1[i][j-1])/(2.*h);
            // Laplacian
            lapV = (v_n[i+1][j]+v_n[i-1][j]+v_n[i][j+1]+v_n[i][j-1]-4.*v_n[i][j])/(h*h);
            // Pressure term
            dpdy = (P[i][j+1]-P[i][j])/h;

            v_star[i][j] = (v_n[i][j] + dt*(-1.5*H_V + 0.5*H_V_old - dpdy + nu*lapV) + (dt/dtau)*ramp*I_V[i][j]*v_s[i][j])/(1.+ramp*I_V[i][j]*dt/dtau);

            /* the value of v_star on the boundaries (j=0, j=n-2) is set to zero at allocation */

        }
    }

#ifdef CAVITY
    for (i=1; i<m-1; i++){
        v_star[i][0] = v_n[i][1]; //(4.*v_n[i][1] - v_n[i][2])/3.;
        v_star[i][n-2] = v_n[i][n-3]; //(4.*v_n[i][n-3] - v_n[i][n-4])/3.;
    }
#endif
}


PetscErrorCode poisson_solver(Data* data, int myrank, int nbproc)
{
    double** u_star = data->u_star;
    double** v_star = data->v_star;
    double** phi = data->phi;

    int M = data->M;
    int N = data->N;
    int m = data->m;
    int n = data->n;
    double h = data->h;
    double dt = data->dt;

    /* Solve the linear system Ax = b for a 2-D poisson equation on a structured grid */
    KSP sles;
    Mat A;
    Vec b, x;

    double div_u_star;
    int r, rowStart, rowEnd, i, j, ii, jj, its;
    int mytag = 12;
    int my_rowStart, my_rowEnd;
    double*  my_array;
    double* array = malloc(M*N*sizeof(double));
    MPI_Status status[3];


    /* Create the Laplacian matrix : A  */
    MatCreate( PETSC_COMM_WORLD, &A );
    MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, M*N, M*N);
    MatSetType(A, MATAIJ);
    MatSeqAIJSetPreallocation(A, 5, NULL); 
    MatMPIAIJSetPreallocation(A, 5, NULL, 5, NULL);
    //PetscInt nz = 5;
    //MatSeqAIJSetPreallocation(A, nz, NULL);
    //MatSetUp(A);
    //MatSetFromOptions(A);
    //MatSetOption(A, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_FALSE);
    MatGetOwnershipRange(A, &rowStart, &rowEnd);
    PetscPrintf(PETSC_COMM_WORLD, "End row is %d \n", rowEnd);

    for(r = rowStart; r<rowEnd; r++){
        ii = r/N; jj=r%N;
        if(ii>0){
            MatSetValue(A, r, r-N, -1., INSERT_VALUES);
        }
        if(jj>0){
            MatSetValue(A, r, r-1, -1., INSERT_VALUES);
        }
        MatSetValue(A, r, r, 4., INSERT_VALUES);
        if(jj<N-1){
            MatSetValue(A, r, r+1, -1., INSERT_VALUES);
        }
        if(ii<M-1){
            MatSetValue(A, r, r+N, -1., INSERT_VALUES);
        }
#ifdef CHANNEL
        if(ii == 0 || jj == 0 || jj == N-1)
        {
            MatSetValue(A, r, r, 3., INSERT_VALUES);
        }
        if(ii == 0 && (jj==0 || jj == N-1) )
        {
            MatSetValue(A, r, r, 2., INSERT_VALUES);
        }

        if(ii == M-1)
        {
            MatSetValue(A, r, r, 5., INSERT_VALUES);
        }
        if(ii == M-1 && (jj==0 || jj==N-1))
        {
            MatSetValue(A, r, r, 4., INSERT_VALUES);
        }
#endif
#ifdef CAVITY
        if(ii == 0 || ii== M-1 || jj == 0 || jj == N-1)
        {
            MatSetValue(A, r, r, 3., INSERT_VALUES);
        }
        if(ii == 0  && (jj==0 || jj == N-1) )
        {
            MatSetValue(A, r, r, 2., INSERT_VALUES);
        }
        if(ii == M-1 && (jj==0 || jj==N-1))
        {
            MatSetValue(A, r, r, 2., INSERT_VALUES);
        }
#endif
    }
    PetscErrorCode  ierr;
    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);


    /* Create the right-hand-side vector : b */
    VecCreate (PETSC_COMM_WORLD, &b );
    VecSetSizes(b, PETSC_DECIDE, M*N);
    VecSetFromOptions(b);
    VecGetOwnershipRange( b, &rowStart, &rowEnd );
    for(r = rowStart; r< rowEnd; r++){
        ii = r/N; jj=r%N;
        i = ii+1; j = jj+1;
        div_u_star = (u_star[i][j]-u_star[i-1][j]+v_star[i][j]-v_star[i][j-1])/h;

        VecSetValue(b, r, -(h*h/dt)*div_u_star, INSERT_VALUES);
    }
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);


    /*Solve the linear system of equations */
    VecDuplicate( b, &x );
    KSPCreate(PETSC_COMM_WORLD, &sles );
    KSPSetOperators(sles, A, A);
    KSPSetFromOptions( sles );
    PetscPrintf(PETSC_COMM_WORLD,"Assembly is done \n");
    KSPSolve( sles, b, x );
    KSPGetIterationNumber( sles, &its );
    PetscPrintf( PETSC_COMM_WORLD, "Solution to Poisson eqn in %d iterations \n ", its);

    /*Transfer the solution to phi[i][j] */
    VecGetArray(x, &my_array);

    if ((M*N) % nbproc == 0){
        int length = ((M*N)/nbproc);
        MPI_Allgather(my_array, length, MPI_DOUBLE, array, length, MPI_DOUBLE, PETSC_COMM_WORLD);
        for (r = 0; r<M*N; r++){
            i = r/N; j = r%N;
            ii = i+1; jj = j+1;
            phi[ii][jj] = array[r];
        }
        // update ghost points on phi
        for(ii=0; ii<m; ii++){
            /* cancel gradients : dp/dn=0 --> dphi/dn = 0*/
            phi[ii][0] = phi[ii][1];
            phi[ii][n-1] = phi[ii][n-2];
        }
        for(jj=0; jj<n; jj++){
            /*inflow : continuity of pressure gradient  */
            phi[0][jj] = phi[1][jj];
#ifdef CHANNEL
            /*outflow : zero pressure at outlet */
            phi[m-1][jj] = -phi[m-2][jj];
#endif
#ifdef CAVITY
            phi[m-1][jj] = phi[m-2][jj];
#endif
        }
    }
    else{
        if (myrank == 0){
            for (r=rowStart; r<rowEnd; r++) {
                array[r] = my_array[r];
            }
            for (int k = 1; k < nbproc; k++){
                MPI_Recv(&my_rowStart, 1, MPI_INT, k, mytag+1, PETSC_COMM_WORLD, &status[1] );
                MPI_Recv(&my_rowEnd, 1, MPI_INT, k, mytag+2, PETSC_COMM_WORLD, &status[2] );
                int length_proc = my_rowEnd - my_rowStart;
                MPI_Recv(my_array, length_proc, MPI_DOUBLE, k, mytag, PETSC_COMM_WORLD, &status[0] ) ;
                int R;
                for (r=0; r<length_proc; r++){
                    R = r + my_rowStart;
                    array[R] = my_array[r];
                }
            }
        }
        else{
            MPI_Send(&rowStart, 1, MPI_INT, 0, mytag+1, PETSC_COMM_WORLD);
            MPI_Send(&rowEnd, 1, MPI_INT, 0, mytag+2, PETSC_COMM_WORLD);
            int length = rowEnd-rowStart;
            MPI_Send(my_array, length, MPI_DOUBLE, 0, mytag, PETSC_COMM_WORLD);
        }
        MPI_Bcast(array, M*N, MPI_DOUBLE, 0, PETSC_COMM_WORLD);

        for(r = 0; r< M*N; r++){
            i = r/N; j = r%N;
            ii = i+1; jj = j+1;
            phi[ii][jj] = array[r];
        }
        // update ghost points on phi
        for(ii=0; ii<m; ii++){
            /* cancel gradients : dp/dn=0 --> dphi/dn = 0*/
            phi[ii][0] = phi[ii][1];
            phi[ii][n-1] = phi[ii][n-2];
        }
        for(jj=0; jj<n; jj++){
            /*inflow : continuity of pressure gradient  */
            phi[0][jj] = phi[1][jj];
#ifdef CHANNEL
            /*outflow : zero pressure at outlet */
            phi[m-1][jj] = -phi[m-2][jj];
#endif
#ifdef CAVITY
            phi[m-1][jj] = phi[m-2][jj];
#endif
        }
    }
    VecRestoreArray(x, &my_array);
    //free(my_array);
    free(array);

    MatDestroy( &A );
    VecDestroy( &b ); VecDestroy( &x );
    KSPDestroy( &sles );

    return ierr;
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
    //double y_ch, u_poiseuille;
    /* Set boundary for u_new, v_new */

    for (i = 0; i < m - 1; i++) {
        for (j = 0; j < n - 1; j++) {
            u_new[i][j] = u_star[i][j] - dt * (phi[i + 1][j] - phi[i][j]) / h;
            v_new[i][j] = v_star[i][j] - dt * (phi[i][j + 1] - phi[i][j]) / h;
            P[i][j] += phi[i][j];
        }
    }

//#ifdef CHANNEL
// INFLOW/OUTFLOW

//    for (i = 1; i < m - 1 ; i++) {
//        //v = 0 everywhere
//        u_new[i][n-2] = u_star[i][n-2] - dt * (phi[i + 1][n-2] - phi[i][n-2]) / h;
//        P[i][n-2] += phi[i][n-2];
//    }
//    for (j = 1; j < n - 1; j++) {
//        u_new[0][j] = u_n[0][j]; //inflow
//        u_new[m-2][j] = u_star[m-2][j] - dt * (phi[m-1][j] - phi[m-2][j]) / h;
//        v_new[m-2][j] = v_star[m-2][j] - dt * (phi[m-2][j + 1] - phi[m-2][j]) / h;
//        P[m-2][j] += phi[m-2][j];
//    }
//#endif

//#ifdef CAVITY
//
//    for (i = 1; i < m - 1 ; i++) {
//        v_new[i][0] = v_star[i][0];
//        v_new[i][n-2] = v_star[i][n-2];
//
//        u_new[i][n-2] = u_star[i][n-2] - dt * (phi[i + 1][n-2] - phi[i][n-2]) / h;
//        P[i][n-2] += phi[i][n-2];
//    }
//    for (j = 1; j < n-1; j++){
//        //SLIP
//        u_new[0][j] = u_star[0][j];
//        u_new[m-2][j] = u_star[m-2][j];
//
//        v_new[m-2][j] = v_star[m-2][j] - dt * (phi[m-2][j + 1] - phi[m-2][j]) / h;
//        P[m-2][j] += phi[m-2][j];
//    }
//
//
//#endif

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

    if (data->ramp > 0) { //data->iter > 1
        free2Darray(T_n_1, m);
        free3Darray(C_n_1, Ns, m);
    }

    data->T_n_1 = T_n;
    data->T_n = T_new;

    data->C_n_1 = C_n;
    data->C_n = C_new;

#endif

    if (data->ramp > 0) { //data->iter > 1
        free2Darray(u_n_1, m);
        free2Darray(v_n_1, m);
    }

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
    double* Sp = data->Sp;
    double* II = data->II;
    double** F = data->F;
    double** G = data->G;
    double** M = data->Mz;
    double g = data->g;
    double** Up = data->Up;
    double** Vp = data->Vp;
    double** Omega_p = data->Omega_p;
    double dt = data->dt;

#ifdef TWO_WAY
    dudt[k] = (23.*F[k][2]-16.*F[k][1]+5.*F[k][0])/(12.*Sp[k]*(rho_r - 1.)) - g;
    Up[k][3] = Up[k][2] + dt*dudt[k];
    dvdt[k] = (23.*G[k][2]-16.*G[k][1]+5.*G[k][0])/(12.*Sp[k]*(rho_r - 1.));
    Vp[k][3] = Vp[k][2] + dt*dvdt[k];
    domegadt[k] = (23.*M[k][2]-16.*M[k][1]+5.*M[k][0])/(12.*II[k]*Sp[k]*(rho_r - 1.));
    Omega_p[k][3] = Omega_p[k][2] + dt*domegadt[k];
#endif

    printf("Up = %f \n", Up[k][3]);
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

    compute_Qr(Qr, PP[0][k][2], dH, k);

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
