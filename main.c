#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#include <petsc.h>
#include <petscksp.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscsys.h>
#include "main.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

//#define RECOVER
#define MOVE
//#define TEMP
#define WRITE
#define DISK
//#define ELLIPSE


int main(int argc, char *argv[]){
    
    PetscInitialize(&argc, &argv, 0, 0);
    int rank, nbproc;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank );
    MPI_Comm_size(PETSC_COMM_WORLD, &nbproc);
    printf("Hello from rank %d \n", rank);
    
    
    /* DIMENSIONS */
    Dp = 1;
    ratio_L_d = 3;
    ratio_d_Dp = 3;
    ratio_Dp_h = 30;
    d = ratio_d_Dp*Dp;
    H = d/2.; /* half -height of the channel */
    L = ratio_L_d*d; /* length of the channel is 'ratio' times its width */
    
    /* PHYSICAL PARAMETERS */
    Pr = 0.7;
    Le = 1; /* Lewis number, ratio between Sc and Prandtl */ 
    Sc = Le*Pr; 
    rho_f = 1.;
    rho_p = 1000.;
    rho_r = rho_p/rho_f;
    cp = 1000.;
    cf = 1000.;
    cr = cp/cf;
    alpha_f = nu/Pr;
    Df = make1DDoubleArray(Ns);
    Df[0] = nu/Sc;
    Df[1] = nu/Sc;

    /* FLOW */
    Rp = 100.;
    Um = 1.; 
    Umax = 1.5*Um;  
    nu = Um*Dp/Rp;
    Tm0 = 100. + 273.15; // cup-mixing temperature at the inlet
    Tp0 = 20. + 273.15;
    
    Ns = 2; /* Number of species */
    Np = 1; /* Number of particles */
    
    /* SPECIES PARAMETERS */
    double CA0 = 1.;
    double CB0 = 0.;
    dH = 0;
    
    /* GRID */
    N = ratio_Dp_h*ratio_d_Dp;
    M = ratio_L_d*N;
    n = N + 2; /*for ghost points */
    m = M + 2; /*for ghost points */
    h = d/N;
    
    /* TIME INTEGRATION */
    double CFL = .1; /*Courant-Freidrichs-Lewy condition on convective term */
    double r = .25; /* Fourier condition on diffusive term */
    double dt_CFL = CFL*h/Umax;
    //double dt = CFL*h/Um; 
    double dt_diff = r*h*h/nu;
    dt = fmin(dt_CFL, dt_diff);
    double ratio_dtau_dt = 1e-2;
    dtau = ratio_dtau_dt*dt;
    
    if(rank == 0){
        printf("Rep = %f\n", Rp);
        printf("ratio L/d = %d \n", ratio_L_d);
        printf("ratio d/dp = %d \n", ratio_d_Dp);
        printf("Um = %f\n", Um);
        printf("dt_CFL = %f\n", dt_CFL);
        printf("dt_diff = %f\n", dt_diff);
        //printf("CFL condition set to %f h/Um \n", CFL);
	printf("CFL condition set to %f h/Umax \n", CFL);
        printf("dt = %f\n", dt);
        printf("dtau = %f\n", dtau);
        printf("ratio dtau/dt = %f \n", ratio_dtau_dt);
    }
    int T_write;
    int N_write;
    if (argc >= 3){
        sscanf(argv[1],"%d",&T_write);
        sscanf(argv[2],"%d",&N_write);
    }
    else{
        T_write = 20; /* number of time steps between two writings */
        N_write = 200; /* number of times we write in files */
    }
    double Tf = N_write*T_write*dt;
    t_move = Tf/10.;
    int nKmax = 2;
    int Kmax = 100; /* number of ramping steps */
    
    double t = 0.;
    double t_start;
    int iter, iter_start;
    
    double ramp = 1./Kmax;
    if(rank == 0){
        printf("Write every %d * dt \n", T_write);
        printf("Write %d times \n", N_write);
        printf("Final time : %f \n \n", Tf);
        printf("ramp = %f \n", ramp);
    }
    
    /* PARTICLES POSITION */
    double* xg = calloc(Np,sizeof(double));
    double* yg = calloc(Np,sizeof(double));
    double* dp = calloc(Np,sizeof(double));
    double* rp = calloc(Np,sizeof(double));
    double* theta = calloc(Np,sizeof(double));
    double* Sp = calloc(Np,sizeof(double));
    double surf = 0.; /* particle surface as computed by integration of the mask funciton */
    double* II = calloc(Np,sizeof(double));
    
    xg[0]=H;
    yg[0]=H;
    dp[0]=Dp;
    rp[0]=dp[0]/2.;
    theta[0] = 0; // M_PI/10.
    a = 2.*Dp;
    b = Dp;

    for(int k=0; k<Np; k++){
#ifdef DISK
        Sp[k]=M_PI*rp[k]*rp[k];
        II[k]=(dp[k]*dp[k])/8.; /* IN 2-D !! */
#endif
#ifdef ELLIPSE	    
        Sp[k]=M_PI*a*b;
        II[k]=.25*(a*a+b*b);
#endif
    }
    
#ifdef WRITE
    DEFINE_FILES
    if(rank==0)
    {
        /* Open files to write the results */
        OPEN_FILES
        /* Write DATA */
        WRITE_DATA
    }
#endif
    
    /* Allocate memory */
    MEMORY_ALLOCATION
    
    /* Try to recover previous state */
    OPEN_STATE
    
#ifndef RECOVER
    state_file = NULL;
#endif
    if(state_file){
        RECOVER_STATE
    }
    else{
        
        /**------------------------------- INITIALIZATION of the domain -------------------------------**/
        
        /* VELOCITY : horizontal flow Um  */
        for(int i=0; i<m; i++){
            for(int j=0; j<n; j++){
		double y_ch = (j-0.5)*h - H;
		U[i][j] = Umax*(1.-(y_ch/H)*(y_ch/H));
#ifdef SLIP
                U[i][j] = Um; /*Plug-Flow*/
#endif
                Ustar[i][j] = U[i][j];

                Told[i][j] = Tm0;
            }
        }
        
        /*Initalization of particles temperatures */
        for(int k=0; k<Np; k++){
            Tp[k] = Tp0;
        }
        
        /**------------------------- FIRST STEP WITH EULER EXPLICIT SCHEME ---------------------------**/
        //dt = 5.*dt_CFL; // When static we increase the time step
        
        /* Location PARTICLES */
        
        /* GET MASK FUNCTIONS */
        get_masks(Ip_S, Ip_U, Ip_V, I_S, I_U, I_V, xg, yg, rp, theta, coloring);
        
        /* Initialize solid temperature and field based on initial particle temperature/species concentration */
#ifdef TEMP
        get_Ts(Ts, Ip_S, Tp);
#endif
        /* Get Ghost points */
        get_ghosts(U, V, P, Told, Cold, 0, 0);
        
        /* Solve the flow */
        
        get_Ustar_EE(U, Ustar, V, Aold, Us, P, I_U, ramp);
        get_Vstar_EE(V, Vstar, U, Bold, Vs, P, I_V, ramp);
        
        //old_poisson_solver(Ustar, Vstar, phi, R);
        poisson_solver(Ustar, Vstar, phi, rank, nbproc);
        update_flow(U, V, P, Ustar, Vstar, phi);
        
        /* Solve temperature and specices */
#ifdef TEMP
        update_temp_species_EE(U, V, T, Told, Ts, C_Told, C, Cold, Cs, C_Cold, I_S, ramp);
#endif
        
        /** ------------------------------- RAMPING ------------------------------- **/
        
        /* Starting from a uniform flow without the particles (i.e. ramp = 0 everywhere) and
         then progressively increasing ramp until ramp*chi = 1 at the position of the ellipse.
         At the end of the loop we obtain the steady solution of the flow for a given
         fixed position of the ellipse. */
        
        for(int K=2; K<=nKmax*Kmax; K++){
            /* role de nKmax : atteindre une solution stable avec les particules fixées avant de les "lacher" */
            ramp = fmin(1.,(double) K/Kmax);
            
            PetscPrintf(PETSC_COMM_WORLD, "\n \n ramp = %f\n",ramp);
            
            for (int k=0; k<Np; k++){
                integrate_penalization(F, G, M, QQ, PP, Ip_U, Ip_V, Ip_S, U, V, T, C, Us, Vs, Ts, Cs, xg, yg, rp, Sp, II, k, &surf);
                /* dudt, dvdt, etc = 0 because the particle is fixed */
                compute_forces(Fx, Fy, Tz, F, G, M,  dudt, dvdt, dwdt, Sp, II, k);
                compute_fluxes(Q, Phi, QQ, PP, dTdt, dCdt, Sp, k);
            }
            
            get_ghosts(U, V, P, Told, Cold, 0, 0);
            
            get_Ustar(U, Ustar, V, A, Aold, Us, P, I_U, ramp);
            get_Vstar(V, Vstar, U, B, Bold, Vs, P, I_V, ramp);
            
            poisson_solver(Ustar, Vstar, phi, rank, nbproc);
            //old_poisson_solver(Ustar, Vstar, phi, R);
            update_flow(U, V, P, Ustar, Vstar, phi);
#ifdef TEMP
            update_temp_species(U, V, T, Told, Ts, C_T, C_Told, C, Cold, Cs, C_C, C_Cold, I_S, ramp);
#endif
            
            
        }
        
#ifdef WRITE
        /*INITIAL SOLUTION AFTER RAMPING */
        if(rank==0){
            fprintf(fichier_position, "%3.13e \t %3.13e \t %3.13e \n ", xg[0], yg[0], theta[0]);
            fflush(fichier_position);
            fprintf(fichier_forces, "%3.13e \t %3.13e \t %3.13e \n ", Fx[0], Fy[0], Tz[0]);
            fflush(fichier_forces);
            fprintf(fichier_fluxes, "%3.13e \t %3.13e \n ", Q[0], Phi[0][0]);
            fflush(fichier_fluxes);
            writeFile(fichier_U,U,0,m-1,1,n-1);
            writeFile(fichier_V,V,1,m-1,0,n-1);
            writeFile(fichier_P,P,1,m-1,1,n-1);
            writeFile(fichier_T,T,1,m-1,1,n-1);
            writeFile(fichier_CA,C[0],1,m-1,1,n-1);
            writeFile(fichier_CB,C[1],1,m-1,1,n-1);
            writeFile(fichier_mask,coloring,1,m-1,1,n-1);
            fprintf(fichier_Tp, "%3.13e \n",Tp[0]);
            fflush(fichier_Tp);
        }
#endif
        /** -------------------------------TIME STEPPING FROM BEGINNING ------------------------------- **/
        iter_start = 1;
        t_start = 0.;
    }
    
    
    
    /** -------------------------------TIME STEPPING ------------------------------- **/
    iter = iter_start;
    t = t_start;
    if(rank==0){
        printf(" \n \n iter %d : t = %f \n", iter, t);
    }
    fflush(stdout);
    
    while(t < Tf){
        //if (t >= t_move){
        //    dt = fmin(dt_CFL, dt_diff);
        //}
        /* SOLVE PARTICULAR PHASE */
        int flag_out = 0;
        
        for (int k = 0; k<Np; k++){
            
            /* Integrate penalization term */
            flag_out += integrate_penalization(F, G, M, QQ, PP, Ip_U, Ip_V, Ip_S, U, V, T, C, Us, Vs, Ts, Cs, xg, yg, rp, Sp, II, k, &surf);
            
            /* Velocity - Forces */
            if(t > t_move){
                update_Xp(xg, yg, theta, Up, Vp, Wp, k);
                update_Up(Up, Vp, Wp, dudt, dvdt, dwdt, F, G, M, II, Sp, k);
            }
            compute_forces(Fx, Fy, Tz, F, G, M,  dudt, dvdt, dwdt, Sp, II, k);
            
            /*Temperature - Species - Fluxes */
	    compute_Qr(Qr, PP, k); 
	    update_Tp(Tp, dTdt, QQ, Qr, Sp, k);
            update_Cp(Cp, dCdt, PP, k);
            compute_fluxes(Q, Phi, QQ, PP, dTdt, dCdt, Sp, k);
        }
        
        if (flag_out > 0){
            break;
        }
        /* SOLVE FLUID PHASE */
        
        /*Compute the mask functions */
        get_masks(Ip_S, Ip_U, Ip_V, I_S, I_U, I_V, xg, yg, rp, theta, coloring);
        
        /* Deduce solid velocity field */
        get_Us_Vs(Us, Vs, Ip_U, Ip_V, Up, Vp, Wp, xg, yg);
        
#ifdef TEMP
        get_Ts(Ts, Ip_S, Tp);
        get_Cs(Cs, Ip_S, Cp);
#endif
        get_ghosts(U, V, P, Told, Cold, CA0, CB0);
        
        get_Ustar(U, Ustar, V, A, Aold, Us, P, I_U, ramp);
        get_Vstar(V, Vstar, U, B, Bold, Vs, P, I_V, ramp);
        
        clock_t t_init = clock();
        //old_poisson_solver(Ustar, Vstar, phi, R);
        poisson_solver(Ustar, Vstar, phi, rank, nbproc);
        clock_t t_final = clock();
        double t_Poisson = ((double) (t_final - t_init))/CLOCKS_PER_SEC;
        PetscPrintf(PETSC_COMM_WORLD, "Poisson solver took %f seconds \n", t_Poisson);
        
        update_flow(U, V, P, Ustar, Vstar, phi);
        /*for(int i=0; i<m; i++){
         for(int j=0; j<n; j++){
         if(I_U[i][j]){
         PetscPrintf(PETSC_COMM_WORLD, "U = %1.6e \t vs.\t  Us = %1.6e\n", U[i][j], Us[i][j]);
         }
         }
         }*/
        
#ifdef TEMP
        update_temp_species(U, V, T, Told, Ts, C_T, C_Told, C, Cold, Cs, C_C, C_Cold, I_S, ramp);
#endif
        
        /* Increment time */
        t += dt;
        iter++;
        
        PetscPrintf(PETSC_COMM_WORLD,"\n \n iter %d : t = %f\n", iter, t);
        fflush(stdout);
#ifdef WRITE
        if(rank == 0){
            fprintf(fichier_position, "%3.13e \t %3.13e \t %3.13e \n ", xg[0], yg[0], theta[0]);
            fflush(fichier_position);
            fprintf(fichier_forces, "%3.13e \t %3.13e \t %3.13e \n ", Fx[0], Fy[0], Tz[0]);
            fflush(fichier_forces);
            fprintf(fichier_fluxes, "%3.13e \t %3.13e \n ", Q[0], Phi[0][0]);
            fflush(fichier_fluxes);
	    printf("Surface from outside : %f \n", surf); 
            fprintf(fichier_surface, "%f \n", surf);
            fflush(fichier_surface);
            fprintf(fichier_Tp, "%3.13e \n",Tp[0]);
            fflush(fichier_Tp);
            
            if(iter % T_write == 0){
                writeFile(fichier_U,U,0,m-1,1,n-1);
                writeFile(fichier_V,V,1,m-1,0,n-1);
                writeFile(fichier_P,P,1,m-1,1,n-1);
#ifdef TEMP
                writeFile(fichier_T,T,1,m-1,1,n-1);
                writeFile(fichier_CA,C[0],1,m-1,1,n-1);
                writeFile(fichier_CB,C[1],1,m-1,1,n-1);
#endif
#ifdef MOVE
                writeFile(fichier_mask,coloring,1,m-1,1,n-1);
#endif
                /*Save current state */
                SAVE_STATE
                PetscPrintf(PETSC_COMM_WORLD, "\n State is saved \n");
            }
        }
#endif
        
    }
    printf("rank %d : job done", rank);
    
#ifdef WRITE
    if (rank ==0){
        /* Close files */
        CLOSE_FILES
    }
#endif
    
    /* Free memory */
    FREE_MEMORY
    
    PetscFinalize();
    return 0;
}

void compute_forces(double* Fx, double* Fy, double* Tz, double** F, double** G, double** M, double* dudt, double* dvdt, double* dwdt, double* Sp, double* II,  int k){
    
    //Fx[k] = rho_p*Sp[k]*dudt[k]; 
    //Fy[k] = rho_p*Sp[k]*dvdt[k]; 	
    //Tz[k] = rho_p*II[k]*Sp[k]*dwdt[k];
    PetscPrintf(PETSC_COMM_WORLD,"\n F integration = %1.6e \n", k+1, F[k][2]);
    PetscPrintf(PETSC_COMM_WORLD,"G integration = %1.6e \n", k+1, G[k][2]);
    PetscPrintf(PETSC_COMM_WORLD,"M integration = %1.6e \n", k+1, M[k][2]);

    PetscPrintf(PETSC_COMM_WORLD,"dudt = %1.6e \n", k+1, dudt[k]);
    PetscPrintf(PETSC_COMM_WORLD,"dvdt = %1.6e \n", k+1, dvdt[k]);
    PetscPrintf(PETSC_COMM_WORLD,"domegadt = %1.6e \n", k+1, dwdt[k]);

    Fx[k] = rho_f*(Sp[k]*dudt[k] + F[k][2]);
    Fy[k] = rho_f*(Sp[k]*dvdt[k] + G[k][2]);
    Tz[k] = rho_f*(Sp[k]*II[k]*dwdt[k] + M[k][2]);

    PetscPrintf(PETSC_COMM_WORLD,"Force along -x dir on particle %d = %1.6e [N/m]  \n", k+1, Fx[k]);
    PetscPrintf(PETSC_COMM_WORLD,"Force along -y dir on particle %d = %1.6e [N/m]  \n", k+1, Fy[k]);
    PetscPrintf(PETSC_COMM_WORLD,"Torque on particle %d = %1.6e [N]  \n", k+1, Tz[k]);
}

void compute_fluxes(double* Q, double** Phi, double** QQ, double*** PP, double* dTdt, double** dCdt, double* Sp, int k){
    Q[k] = rho_f*cf*(Sp[k]*dTdt[k] + QQ[k][2]);
    Phi[k][0] = PP[k][0][2];
    PetscPrintf(PETSC_COMM_WORLD,"Heat flux on particle %d = %1.6e [W/m] \n", k+1, Q[k]);
    PetscPrintf(PETSC_COMM_WORLD,"Molar flux of A on particle %d = %1.6e [mol/(m.s)] \n", k+1, Phi[k][0]);
}

void compute_Qr(double** Qr, double*** PP, int k){
	Qr[k][0] = Qr[k][1];
	Qr[k][1] = Qr[k][2];
	Qr[k][2] = PP[k][0][2]*(-dH);
}

int integrate_penalization(double** F, double** G, double** M, double** QQ, double*** PP, double*** Ip_U, double*** Ip_V, double*** Ip_S, double** U, double** V, double** T, double*** C, double** Us, double** Vs,double** Ts, double*** Cs, double* xg, double* yg, double* rp, double* Sp, double* II, int k, double* pointeur_surf)
{
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
    
    int startX = floor((xg[k]-rp[k])/h);
    PetscPrintf(PETSC_COMM_WORLD, "startX = %d \t", startX);
    int endX = ceil((xg[k]+rp[k])/h);
    PetscPrintf(PETSC_COMM_WORLD,"endX = %d \t", endX);
    
    int startY = floor((yg[k]-rp[k])/h);
    PetscPrintf(PETSC_COMM_WORLD,"startY = %d \t", startY);
    int endY = ceil((yg[k]+rp[k])/h);
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
    
    //for(int i=startX; i<=endX; i++){//
    for(int i = 0; i<m; i++){
        double xV = (i-0.5)*h;
        //for(int j=startY; j<=endY; j++){
        for(int j=0; j<n; j++){
            double yU = (j-0.5)*h;
            double f = -Ip_U[k][i][j]*(U[i][j]-Us[i][j]);
            double g = -Ip_V[k][i][j]*(V[i][j]-Vs[i][j]);
#ifdef TEMP
            double q = -Ip_S[k][i][j]*(T[i][j]-Ts[i][j]);
            double* qm = make1DDoubleArray(Ns);
            for(int s=0; s<Ns; s++){
                qm[s] = -Ip_S[k][i][j]*(C[s][i][j]-Cs[s][i][j]);
            }
#endif
            sint += Ip_U[k][i][j]*h*h;

            Fint += f; /* units : m/s^2 */
            Gint += g; /* units : m/s^2 */
            Mint += ((xV-xg[k])*g-(yU-yg[k])*f);/* units: m^2/s^2 */
#ifdef TEMP
            Qint += q*h*h/dtau; /*units : K/s */
            for(int s=0; s<Ns; s++){
                Qmint[s] += qm[s]*h*h/dtau; /*units : mol/m.s */
            }
#endif
        }
    }
    Fint *= h2/dtau; 
    Gint *= h2/dtau; 
    Mint *= h2/dtau; 
    Qint *= h2/dtau; 
  
    F[k][2] = -Fint;
    G[k][2] = -Gint;
    M[k][2] = -Mint;
    QQ[k][2] = -Qint;
    *pointeur_surf = sint; 

#ifdef TEMP
    for(int s=0; s<Ns; s++){
        PP[k][s][2] = -Qmint[k][s][2]*h2/dtau;
    }
    //double Qr = Phi[k][0][2]*Sp[k]*(-dH); /*units : J/m.s */
    //Q[k][2] = -Qint/(Sp[k]*(rho_r*cr-1.)) + Qr/(Sp[k]*(rho_p*cp-rho_f*cf));
#endif
    PetscPrintf(PETSC_COMM_WORLD, "Particle surface is %f\n", *pointeur_surf);
    /*
     // Compute Hydrodynamic forces and fluxes
     
     #ifdef MOVE
     if(t > t_move){
     F_drag[k] = rho_f*(Sp[k]*F[k][2] - Fint);
     //F_drag[k] = rho_f*(-rho_r/(rho_r - 1.))*Fint; //[N/m]
     F_lift[k] = rho_f*(Sp[k]*G[k][2] - Gint);
     //F_lift[k] = rho_f*(-rho_r/(rho_r - 1.))*Gint; //[N/m]
     Torque[k] = rho_f*(Sp[k]*M[2] - Mint); //[N]
     //Torque[k] = rho_f*(-rho_r/(rho_r - 1.))*Mint; //[N]
     #ifdef TEMP
     Q_heat[k] = rho_f*cf*((-rho_r*cr/(rho_r*cr-1.))*Qint + Qr/(rho_p*cp-rho_f*cf)); //[M/m]
     Phi_species[k] = -Qmint[0]; //[mol/s]
     #endif
     }
     else{
     F_drag[k] = -rho_f*Fint; // [N]
     F_lift[k] = -rho_f*Gint;
     Torque[k] = -rho_f*Mint; //[N.m]
     #ifdef TEMP
     Q_heat[k] = rho_f*cf*(-Qint + Qr/(rho_p*cp-rho_f*cf)); // [W]
     Phi_species[k] = -Qmint[0]; //[mol/s]
     #endif
     }
     #endif
     
     #ifndef MOVE
     F_drag[k] = -rho_f*Fint; // [N]
     F_lift[k] = -rho_f*Gint;
     Torque[k] = -rho_f*Mint; //[N.m]
     #ifdef TEMP
     Q_heat[k] = -rho_f*cf*(Qint + Qr/(rho_p*cp-rho_f*cf)); // [W]
     Phi_species[k] = -Qmint[0]; //[mol/s]
     #endif
     // Here, Phi_species corresponds to the flux of A (reactant) at the surface of particle 0 //
     #endif
     PetscPrintf(PETSC_COMM_WORLD,"Force along -x dir on particle %d = %1.6e [N/m]  \n", k+1, F_drag[k]);
     PetscPrintf(PETSC_COMM_WORLD,"Force along -y dir on particle %d = %1.6e [N/m]  \n", k+1, F_lift[k]);
     PetscPrintf(PETSC_COMM_WORLD,"Torque on particle %d = %1.6e [N]  \n", k+1, Torque[k]);
     PetscPrintf(PETSC_COMM_WORLD,"Heat flux on particle %d = %1.6e [W/m] \n", k+1, Q_heat[k]);
     */
    return 0;
}

void get_ghosts(double** U, double** V, double** P, double** Told, double*** Cold, double CA0, double CB0)
{
    /*Ghost points to impose BC's */
    for (int i=0; i<m; i++){
        /* On U */
        /* Bottom Wall : slip : du/dn= 0 */
#ifdef SLIP 
        U[i][0] = U[i][1];
#endif 
	/*No-slip (channel) */
         U[i][0] = -0.2*(U[i][3] - 5.*U[i][2] + 15.*U[i][1]); 
	/* Top Wall : slip : du/dn = 0  */
#ifdef SLIP 
        U[i][n-1] = U[i][n-2];
#endif 
        U[i][n-1] = -0.2*(U[i][n-4] - 5.*U[i][n-3] + 15.*U[i][n-2]);
        
        
#ifdef TEMP
        /* Walls : adiabatic: dTdn = 0, no mass flux */
        Told[i][0] = Told[i][1];
        Told[i][n-1] = Told[i][n-2];
        
        for (int s=0; s<Ns; s++){
            Cold[s][i][0] = Cold[s][i][1];
            Cold[s][i][n-1] = Cold[s][i][n-2];
        }
#endif
    }
    for (int j = 0; j<n; j++){
        /* On V */
        /* Inflow : horizontal flow --> V = 0 */
        V[0][j] = -0.2*(V[3][j] - 5.*V[2][j] + 15.*V[1][j]);
        /* Natural Ouflow : dV/dx =0 */
        V[m-1][j] = V[m-2][j];
        
        /* On P */
        /* Outflow : P = 0 */
        //P[m-1][j] = -P[m-2][j];
#ifdef TEMP
        /* On T and C */
        /* Inflow : T uniform  */
        Told[0][j] = -0.2*(Told[3][j]-5.*Told[2][j]+15.*Told[1][j]-16.*Tm0);
        
        /* Inflow : CA = CA0; CB = CB0 */
        Cold[0][0][j] = -0.2*(Cold[0][3][j]-5.*Cold[0][2][j]+15.*Cold[0][1][j]-16.*CA0);
        Cold[1][0][j] = -0.2*(Cold[1][3][j]-5.*Cold[1][2][j]+15.*Cold[1][1][j]-16.*CB0);
        
        /*Outflow : We cancel axial dispersion d2T/dx2 = 0; d2C/dx2 = 0; */
        Told[m-1][j] = (7.*Told[m-2][j]-5.*Told[m-3][j]+Told[m-4][j])/3.;
        for(int s=0; s<Ns; s++){
            Cold[s][m-1][j] = (7.*Cold[s][m-2][j]-5.*Cold[s][m-3][j]+Cold[s][m-4][j])/3.;
        }
#endif
    }
    PetscPrintf(PETSC_COMM_WORLD, "CA ghost = %f [mol/m^3]\n", Cold[0][0][50]);
}

void get_masks(double*** Ip_S, double*** Ip_U, double*** Ip_V, double** I_S, double** I_U, double** I_V,  double* xg, double* yg, double* rp, double* theta, double** coloring)
{
    for(int i=0; i<m; i++){
        double xU = i*h;
        double xV = (i-0.5)*h;
        for(int j=0; j<n; j++){
            double yU = (j-0.5)*h;
            double yV = j*h;
            double yS = yU;
            double xS = xV;
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

                 Ip_S[k][i][j]=((xV-xg[k])*(xV-xg[k])+(yU-yg[k])*(yU-yg[k])<= rp[k]*rp[k]);
                 Ip_U[k][i][j]=((xU-xg[k])*(xU-xg[k])+(yU-yg[k])*(yU-yg[k])<= rp[k]*rp[k]);
                 Ip_V[k][i][j]=((xV-xg[k])*(xV-xg[k])+(yV-yg[k])*(yV-yg[k])<= rp[k]*rp[k]);
                 
		 double xloc = xV-xg[k];
                 double yloc = yU-yg[k];
                 double delta = 0;
                 coloring[i][j] += Ip_S[k][i][j];
                 if (xloc>0 && yloc>0){ //1st Quadrant
                 delta = atan(yloc/xloc);
                 }
                 if (xloc>0 && yloc<0){ //4th Quadrant
                 delta = atan(yloc/xloc) + 2.*M_PI;
                 }
                 if (xloc<0){
                 delta = atan(yloc/xloc) + M_PI;
                 }
                 if (xloc == 0 && yloc > 0){
                 delta = M_PI/2.;
                 }
                 if (xloc == 0 && yloc < 0){
                 delta = 3.*M_PI/2.;
                 }
                 
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

void get_Cs(double*** Cs, double*** Ip_S, double** Cp)
{
    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            Cs[1][i][j] = 0.;
            for(int k=0; k<Np; k++){
                Cs[1][i][j] += Ip_S[k][i][j]*Cp[k][1];
            }
        }
    }
}
void get_Ts(double** Ts, double*** Ip_S, double* Tp){
    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            Ts[i][j] = 0.;
            for(int k = 0; k<Np; k++){
                Ts[i][j]+= Ip_S[k][i][j]*Tp[k];
                
            }
        }
    }
}

void get_Us_Vs(double** Us, double** Vs, double*** Ip_U, double*** Ip_V, double** Up, double** Vp, double** Wp, double* xg, double* yg){
    for(int i=0; i<m; i++){
        double xV = (i-0.5)*h;
        for(int j=0; j<n; j++){
            double yU = (j-0.5)*h;
            Us[i][j] = 0.;
            Vs[i][j] = 0.;
            for (int k = 0; k<Np; k++){
                Us[i][j]+= Ip_U[k][i][j]*(Up[k][3] - Wp[k][3]*(yU-yg[k]));
                Vs[i][j]+= Ip_V[k][i][j]*(Vp[k][3] + Wp[k][3]*(xV-xg[k]));
            }
        }
    }
}


void get_Ustar(double** U, double** Ustar, double** V, double** A, double** Aold, double** Us, double** P, double** I_U, double ramp)
{
    /* Ustar ADAMS-BASHFORTH 2 */
    for (int i=1; i<m-2; i++){
        for (int j=1; j<n-1; j++){
            /*only interior points*/
            /*Ustar */
            double Vij = .25*(V[i][j]+V[i+1][j]+V[i][j-1]+V[i+1][j-1]); /* average of the four neighbour points */
            A[i][j]= U[i][j]*(U[i+1][j]-U[i-1][j])/(2.*h) + Vij*(U[i][j+1]-U[i][j-1])/(2.*h);
            Ustar[i][j] = (U[i][j] + dt*(-1.5*A[i][j]+0.5*Aold[i][j] - (P[i+1][j]-P[i][j])/h
                                         + nu*(U[i+1][j]+U[i-1][j]+U[i][j+1]+U[i][j-1]-4.*U[i][j])/(h*h)
                                         + ramp*I_U[i][j]*Us[i][j]/dtau))/(1.+ramp*I_U[i][j]*dt/dtau);
            
            Aold[i][j]=A[i][j];
        }
    }
    /*Outflow condition */
    for (int j=1; j<n-1; j++){
	    double y_ch = (j-0.5)*h - H; /*y_ch : channel definition of y-coord.  */
	    double Upoiseuille = Umax*(1.-(y_ch/H)*(y_ch/H));
	    Ustar[m-2][j] = U[m-2][j] - dt*Upoiseuille*(U[m-2][j]-U[m-3][j])/h;
    }
    /*Ustar[0][j] (inflow) is fixed once for all at the beginning */
}

void get_Ustar_EE(double** U, double** Ustar, double** V, double** Aold, double** Us, double** P, double** I_U, double ramp)
{
    /* Ustar EULER EXPLICIT */\
    for (int i=1; i<m-2; i++){
        for (int j=1; j<n-1; j++){
            /*only interior points*/
            /*Ustar */
            double Vij = .25*(V[i][j]+V[i+1][j]+V[i][j-1]+V[i+1][j-1]); /* average of the four neighbour points */
            Aold[i][j]= U[i][j]*(U[i+1][j]-U[i-1][j])/(2.*h) + Vij*(U[i][j+1]-U[i][j-1])/(2.*h);
            Ustar[i][j] = (U[i][j] + dt*(-Aold[i][j] - (P[i+1][j]-P[i][j])/h
                                         + nu*(U[i+1][j]+U[i-1][j]+U[i][j+1]+U[i][j-1]-4.*U[i][j])/(h*h)
                                         + ramp*I_U[i][j]*Us[i][j]/dtau))/(1.+ramp*I_U[i][j]*dt/dtau);
        }
    }
    /*Outflow condition */
    for (int j=1; j<n-1; j++){
        Ustar[m-2][j] = U[m-2][j] - dt*Um*(U[m-2][j]-U[m-3][j])/h;
    }
    /*Ustar[0][j] (inflow) is fixed once for all at the beginning --> but it does not matter with the MAC mesh !*/
}

void get_Vstar(double** V, double** Vstar, double** U, double** B, double** Bold, double** Vs, double** P, double** I_V, double ramp)
{
    /* Vstar  ADAMS-BASHFORTH 2 */
    for (int i=1; i<m-1; i++){
        for (int j=1; j<n-2; j++){
            double Uij = .25*(U[i][j]+U[i-1][j]+U[i][j+1]+U[i-1][j+1]); /* average of the four neighbour points */
            B[i][j] = Uij*(V[i+1][j]-V[i-1][j])/(2.*h) + V[i][j]*(V[i][j+1]-V[i][j-1])/(2.*h);
            Vstar[i][j] = (V[i][j] + dt*(-1.5*B[i][j]+0.5*Bold[i][j] - (P[i][j+1]-P[i][j])/h
                                         + nu*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]-4.*V[i][j])/(h*h)
                                         + ramp*I_V[i][j]*Vs[i][j]/dtau))/(1.+ramp*I_V[i][j]*dt/dtau);
            Bold[i][j]=B[i][j];
            /* the value of Vstar on the boundaries (j=0, j=n-2) does not matter with the MAC mesh */
        }
    }
}

void get_Vstar_EE(double** V, double** Vstar, double** U, double** Bold, double** Vs, double** P, double** I_V, double ramp)
{
    /* Vstar EULER EXPLICIT */
    for (int i=1; i<m-1; i++){
        for (int j=1; j<n-2; j++){
            double Uij = .25*(U[i][j]+U[i-1][j]+U[i][j+1]+U[i-1][j+1]); /* average of the four neighbour points */
            Bold[i][j] = Uij*(V[i+1][j]-V[i-1][j])/(2.*h) + V[i][j]*(V[i][j+1]-V[i][j-1])/(2.*h);
            Vstar[i][j] = (V[i][j] + dt*(-Bold[i][j] - (P[i][j+1]-P[i][j])/h
                                         + nu*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]-4.*V[i][j])/(h*h)
                                         + ramp*I_V[i][j]*Vs[i][j]/dtau))/(1.+ramp*I_V[i][j]*dt/dtau);
        }
    }
}

void poisson_solver(double** Ustar, double** Vstar, double **phi, int myrank, int nbproc)
{
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

    /* CAREFUL : this implementation of POISSON solver is only valid with no-slip boundary condition : V = U = 0 at the walls */ 

    MatCreate( PETSC_COMM_WORLD, &A );
    MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, M*N, M*N);
    //MatSetType(A, MATMPIAIJ);
    //MatMPIAIJSetPreallocation(A, 1, NULL, 4, NULL);
    MatSetUp(A);
    //MatSetFromOptions(A);
    //MatSetOption(A, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_FALSE);
    MatGetOwnershipRange(A, &rowStart, &rowEnd);
    //MatAIJSetPreallocation(A, 5, NULL);
    
    //PetscPrintf(PETSC_COMM_WORLD, "starting row is %d \n", rowStart);
    //printf("starting row is %d \n", rowStart);
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
    }
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    
    /* Create the right-hand-side vector : b */
    VecCreate (PETSC_COMM_WORLD, &b );
    VecSetSizes(b, PETSC_DECIDE, M*N);
    VecSetFromOptions(b);
    VecGetOwnershipRange( b, &rowStart, &rowEnd );
    for(r = rowStart; r< rowEnd; r++){
        ii = r/N; jj=r%N;
        i = ii+1; j = jj+1;
        div_u_star = (Ustar[i][j]-Ustar[i-1][j]+Vstar[i][j]-Vstar[i][j-1])/h;
        
        VecSetValue(b, r, -(h*h/dt)*div_u_star, INSERT_VALUES);
    }
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);
    
    /*Solve the linear system of equations */
    VecDuplicate( b, &x );
    KSPCreate(PETSC_COMM_WORLD, &sles );
    KSPSetOperators( sles, A, A );
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
        for(ii=0; ii<m; ii++){
            /* cancel gradients : dp/dn=0 --> dphi/dn = 0*/
            phi[ii][0] = phi[ii][1];
            phi[ii][n-1] = phi[ii][n-2];
        }
        for(jj=0; jj<n; jj++){
            /*inflow : continuity of pressure gradient  */
            phi[0][jj] = phi[1][jj];
            /*outflow : zero pressure at outlet */
            phi[m-1][jj] = -phi[m-2][jj];
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
                int length_k = my_rowEnd - my_rowStart;
                MPI_Recv(my_array, length_k, MPI_DOUBLE, k, mytag, PETSC_COMM_WORLD, &status[0] ) ;
                for (r=0; r<length_k; r++){
                    int R = r+my_rowStart;
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
        for(ii=0; ii<m; ii++){
            /* cancel gradients : dp/dn=0 --> dphi/dn = 0*/
            phi[ii][0] = phi[ii][1];
            phi[ii][n-1] = phi[ii][n-2];
        }
        for(jj=0; jj<n; jj++){
            /*inflow : continuity of pressure gradient  */
            phi[0][jj] = phi[1][jj];
            /*outflow : zero pressure at outlet */
            phi[m-1][jj] = -phi[m-2][jj];
        }
    }
    VecRestoreArray(x, &my_array);
    free(my_array);
    free(array);
    MatDestroy( &A );
    VecDestroy( &b ); VecDestroy( &x );
    KSPDestroy( &sles );
    
}

void old_poisson_solver(double** Ustar, double** Vstar, double **phi, double** R)
{
    double e = 0.;
    int SORiter = 0;
    do{
        for (int i=1; i<m-1; i++){
            for (int j=1; j<n-1; j++){
                double phistar = .25*(-(h/dt)*(Ustar[i][j]-Ustar[i-1][j]+Vstar[i][j]-Vstar[i][j-1])
                                      + phi[i+1][j] + phi[i-1][j] + phi[i][j+1] + phi[i][j-1]);
                phi[i][j] = alpha*phistar + (1.-alpha)*phi[i][j];
            }
        }
        /*Adapt ghost points on Phi */
        for (int i=1; i<m-1; i++){
            phi[i][0] = phi[i][1];
            phi[i][n-1] = phi[i][n-2];
        }
        
        /*LEFT AND RIGHT */
        for (int j=1; j<n-1; j++){
            phi[0][j] = phi[1][j];
            phi[m-1][j] = -phi[m-2][j];
        }
        
        double Sum = 0.;
        
        for (int i=1; i<m-1; i++){
            for (int j=1; j<n-1; j++){
                R[i][j] = -(Ustar[i][j]-Ustar[i-1][j]+Vstar[i][j]-Vstar[i][j-1])/(h*dt)
                +(phi[i+1][j]+phi[i-1][j]+phi[i][j+1]+phi[i][j-1]-4.*phi[i][j])/(h*h);
                /* CACLULER LA SOMME */
                Sum += R[i][j]*R[i][j];
                
            }
        }
        
        e = (h*dt*d/Um)*pow(Sum/(L*d),0.5);
        SORiter ++;
    } while (e > SORtol && SORiter < SORitermax);
    PetscPrintf(PETSC_COMM_WORLD, "e = %1.3e after %d/%d iterations \n",e,SORiter,SORitermax);
    
}

void reset_phi(double** phi)
{
    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            phi[i][j] = 0.;
        }
    }
}


void update_flow(double** U, double** V, double** P, double** Ustar, double** Vstar, double** phi)
{
    for (int i=1; i<m-1; i++){
        int j;
        for (j=1; j<n-2; j++){
            U[i][j] = Ustar[i][j] - dt*(phi[i+1][j]-phi[i][j])/h;
            P[i][j] += phi[i][j];
            V[i][j] = Vstar[i][j] - dt*(phi[i][j+1]-phi[i][j])/h;
        }
        /* Last value of j is n-2 */
        U[i][j] = Ustar[i][j] - dt*(phi[i+1][j]-phi[i][j])/h;
        P[i][j] += phi[i][j];
      
	V[i][0] = 0; 
	V[i][n-2] = 0; 
#ifdef SLIP
       	V[i][0] = (4.*V[i][1]-V[i][2])/3.;/* to cancel gradient at the slip boundary : dv/dn = 0 */  /* Second-order accuracy */
        V[i][n-2] = (4.*V[i][n-3]-V[i][n-4])/3.; /* to cancel gradient at the slip boundary : dv/dn = 0 */
#endif
    }
}

void update_temp_species(double** U, double** V,  double** T, double** Told, double** Ts, double** C_T, double** C_Told, double*** C, double*** Cold, double*** Cs, double*** C_C, double*** C_Cold, double** I_S, double ramp)
{
    /* Adams-Bashfort 2 */
    for (int i=1; i<m-1; i++){
        for (int j=1; j<n-1; j++){
            /*only interior points*/
            /* T */
            double Uij = .5*(U[i][j]+U[i-1][j]); /* average of the two neighbour points along x-direction */
            double Vij = .5*(V[i][j]+V[i-1][j]);
            C_T[i][j]= Uij*(Told[i+1][j]-Told[i-1][j])/(2.*h) + Vij*(Told[i][j+1]-Told[i][j-1])/(2.*h);
            T[i][j] = (Told[i][j] + dt*(-1.5*C_T[i][j] + 0.5*C_Told[i][j]
                                        + alpha_f*(Told[i+1][j]+Told[i-1][j]+Told[i][j+1]+Told[i][j-1]-4.*Told[i][j])/(h*h)
                                        + ramp*I_S[i][j]*Ts[i][j]/dtau))/(1.+ramp*I_S[i][j]*dt/dtau);
            for (int s=0; s<Ns; s++){
                C_C[s][i][j] = Uij*(Cold[s][i+1][j]-Cold[s][i-1][j])/(2.*h) + Vij*(Cold[s][i][j+1]-Cold[s][i][j-1])/(2.*h);
                C[s][i][j] = (Cold[s][i][j] + dt*(-1.5*C_C[s][i][j] + 0.5*C_Cold[s][i][j]
                                                  + Df[s]*(Cold[s][i+1][j]+Cold[s][i-1][j]+Cold[s][i][j+1]+Cold[s][i][j-1]-4.*Cold[s][i][j])/(h*h)
                                                  + ramp*I_S[i][j]*Cs[s][i][j]/dtau))/(1.+ramp*I_S[i][j]*dt/dtau);
            }
        }
    }
    for (int i=1; i<m-1; i++){
        for (int j=1; j<n-1; j++){
            C_Told[i][j] = C_T[i][j];
            Told[i][j] = T[i][j];
            for (int s=0; s<Ns; s++){
                C_Cold[s][i][j] = C_C[s][i][j];
                Cold[s][i][j] = C[s][i][j];
            }
        }
    }
}

void update_temp_species_EE(double** U, double** V, double** T, double** Told, double** Ts, double** C_Told, double*** C, double*** Cold, double*** Cs, double*** C_Cold, double** I_S, double ramp)
{
    /* Euler Explicit */
    for (int i=1; i<m-1; i++){
        for (int j=1; j<n-1; j++){
            /*only interior points*/
            /* T */
            double Uij = .5*(U[i][j]+U[i-1][j]); /* average of the two neighbour points along x-direction */
            double Vij = .5*(V[i][j]+V[i-1][j]);
            C_Told[i][j]= Uij*(Told[i+1][j]-Told[i-1][j])/(2.*h) + Vij*(Told[i][j+1]-Told[i][j-1])/(2.*h);
            T[i][j] = (Told[i][j] + dt*(-C_Told[i][j]
                                        + alpha_f*(Told[i+1][j]+Told[i-1][j]+Told[i][j+1]+Told[i][j-1]-4.*Told[i][j])/(h*h)
                                        + ramp*I_S[i][j]*Ts[i][j]/dtau))/(1.+ramp*I_S[i][j]*dt/dtau);
            for (int s=0; s<Ns; s++){
                C_Cold[s][i][j] = Uij*(Cold[s][i+1][j]-Cold[s][i-1][j])/(2.*h) + Vij*(Cold[s][i][j+1]-Cold[s][i][j-1])/(2.*h);
                C[s][i][j] = (Cold[s][i][j] + dt*(-C_Cold[s][i][j]
                                                  + Df[s]*(Cold[s][i+1][j]+Cold[s][i-1][j]+Cold[s][i][j+1]+Cold[s][i][j-1]-4.*Cold[s][i][j])/(h*h)
                                                  + ramp*I_S[i][j]*Cs[s][i][j]/dtau))/(1.+ramp*I_S[i][j]*dt/dtau);
            }
        }
    }
    for (int i=1; i<m-1; i++){
        for (int j=1; j<n-1; j++){
            Told[i][j] = T[i][j];
            for (int s=0; s<Ns; s++){
                Cold[s][i][j] = C[s][i][j];
            }
        }
    }
}

void update_Xp(double* xg, double* yg, double* theta, double** Up, double** Vp, double** Wp,  int k)
{
    xg[k] += dt*(23.*Up[k][2]-16.*Up[k][1]+5.*Up[k][0])/12.;
    yg[k] += dt*(23.*Vp[k][2]-16.*Vp[k][1]+5.*Vp[k][0])/12.;
    theta[k] += dt*(23.*Wp[k][2]-16.*Wp[k][1]+5.*Wp[k][0])/12.;
    PetscPrintf(PETSC_COMM_WORLD,"Position of the center of mass of particle %d: (x,y) = (%f,%f) \n", k+1, xg[k], yg[k]);
    PetscPrintf(PETSC_COMM_WORLD,"Angle: theta  = %f \n", theta[k]);
    
}

void update_Up(double** Up, double** Vp, double** Wp, double* dudt, double* dvdt, double* domegadt, double** F, double** G, double** M, double* II, double* Sp, int k)
{
#ifdef DISK
    dudt[k] = (23.*F[k][2]-16.*F[k][1]+5.*F[k][0])/(12.*Sp[k]*(rho_r - 1.)) ;
    Up[k][3] = Up[k][2] + dt*dudt[k];
    dvdt[k] = (23.*G[k][2]-16.*G[k][1]+5.*G[k][0])/(12.*Sp[k]*(rho_r - 1.));
    Vp[k][3] = Vp[k][2] + dt*dvdt[k];
#endif 
    domegadt[k] = (23.*M[k][2]-16.*M[k][1]+5.*M[k][0])/(12.*II[k]*Sp[k]*(rho_r - 1.));
    Wp[k][3] = Wp[k][2] + dt*domegadt[k];
    
    Up[k][0] = Up[k][1]; Up[k][1] = Up[k][2]; Up[k][2] = Up[k][3];
    Vp[k][0] = Vp[k][1]; Vp[k][1] = Vp[k][2]; Vp[k][2] = Vp[k][3];
    Wp[k][0] = Wp[k][1]; Wp[k][1] = Wp[k][2]; Wp[k][2] = Wp[k][3];
}


void update_Tp(double* Tp, double* dTdt, double** QQ, double** Qr, double* Sp, int k)
{
    dTdt[k] = (23.*QQ[k][2]-16.*QQ[k][1]+5.*QQ[k][0])/(12.*Sp[k]*(rho_r*cr - 1.)) + (23.*Qr[k][2]-16.*Qr[k][1]+5.*Qr[k][0])/(12.*Sp[k]*(rho_p*cp - rho_f*cf));
    Tp[k] += dt*dTdt[k];
    PetscPrintf(PETSC_COMM_WORLD,"Temperature of particle %d: Tp = %f[K] \n", k+1, Tp[k]);
}

void update_Cp(double** Cp, double** dCdt, double*** PP, int k)
{
    dCdt[k][1] = (23.*(PP[k][0][2]+PP[k][1][2])-16.*(PP[k][0][1]+PP[k][1][1])+5.*(PP[k][0][0]+PP[k][1][0]))/12.;
    Cp[k][0] = 0.;
    Cp[k][1] += dt*dCdt[k][1];
    //printf("Concentration of B in particle %d: Cp_B = %3.13e [mol/m3] \n", k+1, Cp[k][1]);
}

/* I/O related functions */

/*
 void combine(char* destination, const char* path1, const char* path2)
 {
 if(path1 == NULL && path2 == NULL) {
 strcpy(destination, "");;
 }
 else if(path2 == NULL || strlen(path2) == 0) {
 strcpy(destination, path1);
 }
 else if(path1 == NULL || strlen(path1) == 0) {
 strcpy(destination, path2);
 }
 else {
 char directory_separator[] = "/";
 const char *last_char = path1;
 while(*last_char != '\0')
 last_char++;
 int append_directory_separator = 0;
 if(strcmp(last_char, directory_separator) != 0) {
 append_directory_separator = 1;
 }
 strcpy(destination, path1);
 if(append_directory_separator)
 strcat(destination, directory_separator);
 strcat(destination, path2);
 }
 }
 */


void writeFile(FILE* file, double **data, int iStart, int iEnd, int jStart, int jEnd)
{
    if (file != NULL){
        for (int i = iStart; i<iEnd; i++){
            int j;
            for (j = jStart; j<jEnd-1; j++) {
                fprintf(file, "%3.13e \t",data[i][j]);
            }
            fprintf(file, "%3.13e", data[i][j]);
            fprintf(file, "\n");
            fflush(file);
        }
    }
    else{
        printf(" An error occured in writeFile : invalid file pointer.\n");
    }
}

double* make1DDoubleArray(int arraySize) {
    double* theArray = (double*) calloc(arraySize, sizeof(double));
    return theArray;
}

double** make2DDoubleArray(int arraySizeX, int arraySizeY) {
    double** theArray;
    theArray = (double**) malloc(arraySizeX*sizeof(double*));
    for (int ix = 0; ix < arraySizeX; ix++){
        theArray[ix] =(double*) calloc(arraySizeY, sizeof(double));
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
            theArray[k][ix] = calloc(arraySizeY, sizeof(double));
        }
    }
    return theArray;
}

int** make2DIntArray(int arraySizeX, int arraySizeY) {
    int** theArray;
    theArray = (int**) malloc(arraySizeX*sizeof(int*));
    for (int ix = 0; ix < arraySizeX; ix++) {
        theArray[ix] =(int*) calloc(arraySizeY, sizeof(int));
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
            theArray[k][ix] = calloc(arraySizeY, sizeof(int));
        }
    }
    return theArray;
}
