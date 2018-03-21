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

#define WRITE


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

    data.Dp = 1;
    data.L = 2.2;
    data.H = 0.5*data.L;
    data.Omega = 1.25;
    data.R1 = 0.4;
    data.R2 = 1; 

    data.h = data.L/200.;
    data.eps = 4.*data.h;

    /* NON-DIMENSIONAL NUMBERS */
    data.Pr = 0.7;
    data.Le = 1; /* Lewis number, ratio between Sc and Prandtl */
    data.Sc = data.Le*data.Pr;
    data.Rep = 100.;
    data.Fr = sqrt(1e3);

    /* FLOW */
    data.u_m = data.Omega*data.R2;
    data.nu = 0.1; //data.u_m*data.Dp/data.Rep;

    /* ENERGY */
    data.alpha_f = data.nu/data.Pr;
    data.Tm0 = 1; // 100. + 273.15; // cup-mixing temperature at the inlet
    data.Tp0 = 0.5; //20. + 273.15; // initial particle temperature

    /* PHYSICAL PARAMETERS */
    data.rho_f = 1.;
    data.rho_p = 1000.;
    data.rho_r = data.rho_p/data.rho_f;
    data.cp = 1000.;
    data.cf = 1000.;
    data.cr = data.cp/data.cf;

    /* GRID */
    data.M = (int) (data.L/data.h);
    data.N = data.M;
    data.n = data.N + 2; /*for ghost points */
    data.m = data.M + 2; /*for ghost points */
    data.Np = 1;
    data.Ns = 2;

    /* TIME INTEGRATION */
    data.CFL = .1; /*Courant-Freidrichs-Lewy condition on convective term */
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
        data.T_write = 30; /* number of time steps between two writings */
        data.N_write = 100; /* number of times we write in files */
    }


    double Tf = data.N_write*data.T_write*data.dt;
    data.Tf = Tf;
    data.t_move = 0; //data.Tf/10.;
    data.nKmax = 1;
    data.Kmax = 100; /* number of ramping steps */

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

    data.coloring = make2DDoubleArray(m,n);
    data.I_S = make2DDoubleArray(m,n);
    data.I_U = make2DDoubleArray(m,n);
    data.I_V = make2DDoubleArray(m,n);

    data.omega = make2DDoubleArray(m,n);
    data.phi = make2DDoubleArray(m,n);
    data.P = make2DDoubleArray(m,n);
    data.Reh = make2DDoubleArray(m,n);
    data.Reh_omega = make2DDoubleArray(m,n);

    data.u_n = make2DDoubleArray(m,n);
    data.u_n_1 = make2DDoubleArray(m,n);
    data.u_s = make2DDoubleArray(m,n);
    data.u_star = make2DDoubleArray(m,n);

    data.v_n = make2DDoubleArray(m,n);
    data.v_n_1 = make2DDoubleArray(m,n);
    data.v_s = make2DDoubleArray(m,n);
    data.v_star = make2DDoubleArray(m,n);



    /** -------- Some files creation and data writing-------- **/

    //FILE* fichier_data = fopen("results/data.txt", "w");
    //writeData(fichier_data, data);
    //fclose(fichier_data);

    FILE* fichier_stat = fopen("results/stats.txt", "w");
    FILE* fichier_forces = fopen("results/forces.txt", "w");
    FILE* fichier_fluxes = fopen("results/fluxes.txt", "w");
    FILE* fichier_particle = fopen("results/particle.txt", "w");

    /** -------------------------------TIME STEPPING ------------------------------- **/

    get_masks(&data);
    get_Us_Vs(&data);

    for(int K=0; K<=data.nKmax*data.Kmax; K++) {
        data.ramp = fmin(1., (double) K / data.Kmax);
        PetscPrintf(PETSC_COMM_WORLD, "\n \n BEGIN ramp = %f \n", data.ramp);
        get_ghosts(&data, data.Tm0, data.C0);
        get_Ustar_Vstar(&data, data.ramp);
        poisson_solver(&data, rank, nbproc);
        update_flow(&data);
    }

#ifdef WRITE
    /*INITIAL SOLUTION (t=0) AFTER RAMPING */
    if(rank==0){
        writeFields(&data, 0);
    }
#endif

    /** -------------------------------TIME STEPPING FROM BEGINNING ------------------------------- **/
    iter_start = 1;
    t_start = 0.;


    while(t < data.Tf){

        PetscPrintf(PETSC_COMM_WORLD, "\n \n BEGIN iter %d : t = %f \n", data.iter, t);

        /** --- SOLVE FLUID PHASE --- */
        get_ghosts(&data, data.Tm0, data.C0);
        get_Ustar_Vstar(&data, data.ramp);

        clock_t t_init = clock();
        poisson_solver(&data, rank, nbproc);
        clock_t t_final = clock();
        double t_Poisson = ((double) (t_final - t_init))/CLOCKS_PER_SEC;
        PetscPrintf(PETSC_COMM_WORLD, "Poisson solver took %f seconds \n", t_Poisson);

        update_flow(&data);


#ifdef WRITE
        if(rank == 0){
            if(data.iter % data.T_write == 0){
                writeFields(&data, data.iter);
            }
        }
#endif
        /* Increment time */
        t += data.dt;
        data.iter ++;

    }


    /* Free memory */
    free2Darray(data.u_n,m), free2Darray(data.u_n_1,m), free2Darray(data.u_star,m), free2Darray(data.u_s,m);
    free2Darray(data.v_n,m), free2Darray(data.v_n_1,m), free2Darray(data.v_star,m), free2Darray(data.v_s,m);
    free2Darray(data.omega, m); free2Darray(data.Reh,m); free2Darray(data.Reh_omega,m);
    free2Darray(data.P,m), free2Darray(data.phi, m);
    free2Darray(data.I_S, m), free2Darray(data.I_U, m), free2Darray(data.I_V, m), free2Darray(data.coloring, m);

    PetscFinalize();
    return 0;
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


void get_ghosts(Data* data, double T0, double* C0)
{
    double** u_n = data->u_n;
    double** v_n = data->v_n;

    int m = data->m;
    int n = data->n;

    /*Ghost points to impose BC's */
    /* Along y = -H and y = H */
    for (int i=0; i<m; i++){
        /* On u_n */

        /* No-slip (channel) */
        u_n[i][0] = -0.2*(u_n[i][3] - 5.*u_n[i][2] + 15.*u_n[i][1]);
        u_n[i][n-1] = -0.2*(u_n[i][n-4] - 5.*u_n[i][n-3] + 15.*u_n[i][n-2]);

    }
	/* Along x = 0 and x = L */
    for (int j = 0; j<n; j++){
        /* On v_n */
        /* no-slip */
        v_n[0][j] = -0.2*(v_n[3][j] - 5.*v_n[2][j] + 15.*v_n[1][j]);
        /* Natural outflow : dV/dx =0 */
        v_n[m-1][j] = -0.2*(v_n[m-4][j] - 5.*v_n[m-1][j] + 15.*v_n[m-2][j]);
    }
}

void get_masks(Data* data)
{
    double** I_S = data->I_S;
    double** I_U = data->I_U;
    double** I_V = data->I_V;
    double** coloring = data->coloring;

    int m = data->m;
    int n = data->n;
    double h = data->h;
    double H = data->H;

    double xU, xV, yU, yV, yS, xS;
    double xloc, yloc, rloc;

    for(int i=0; i<m; i++){
        xU = i*h;
        xV = (i-0.5)*h;
        for(int j=0; j<n; j++){
            yU = (j-0.5)*h;
            yV = j*h;
            yS = yU;
            xS = xV;

            xloc = xS-H; yloc = yS-H;
            rloc = sqrt(xloc*xloc + yloc*yloc);
            I_S[i][j]=(rloc<data->R1 || rloc>data->R2);

            xloc = xU-H; yloc = yU-H;
            rloc = sqrt(xloc*xloc + yloc*yloc);
            I_U[i][j]=(rloc<data->R1 || rloc>data->R2);

            xloc = xV-H; yloc = yV-H;
            rloc = sqrt(xloc*xloc + yloc*yloc);
            I_V[i][j]=(rloc<data->R1 || rloc>data->R2);

            coloring[i][j] = I_S[i][j];
        }
    }
}

void get_Us_Vs(Data* data){

    double** u_s = data-> u_s;
    double** v_s = data-> v_s;

    int m = data->m;
    int n = data->n;
    double h = data->h;
    double H = data->H;

    double xU, xV, yU, yV; 
    for(int i=0; i<m; i++){
        xV = (i-H)*h;
	    xU = i*h;
        for(int j=0; j<n; j++){
            yU = (j-H)*h;
	        yV = j*h;
            if((xU-H)*(xU-H) + (yU-H)*(yU-H) > data->R2*data->R2)
            {
                u_s[i][j]= -data->Omega*(yU-H);
            }
            if ((xV - H) * (xV - H) + (yV - H) * (yV - H) > data->R2 * data->R2)
            {
                v_s[i][j] = data->Omega*(xV - H);
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

            if (H_U != H_U_old && data->ramp == 0){
                printf("Problem at first time step ! ");
            }


            u_star[i][j] = (u_n[i][j] + dt*(-1.5*H_U + 0.5*H_U_old - dpdx + nu*lapU) + (dt/dtau)*ramp*I_U[i][j]*u_s[i][j])/(1.+ramp*I_U[i][j]*dt/dtau);
        }
    }
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
#ifndef TEMP
    if (data->iter > 1) {
        free2Darray(u_n_1, m);
        free2Darray(v_n_1, m);
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
        if(ii == 0 || jj == 0 || jj == N-1 || ii == M-1)
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
            /*outflow : zero pressure at outlet */
            phi[m-1][jj] = phi[m-2][jj];
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
            /*outflow : zero pressure at outlet */
            phi[m-1][jj] = phi[m-2][jj];
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

    double** u_n = data->u_n;
    double** v_n = data->v_n;
    double** P = data->P;

    double** u_star = data->u_star;
    double** v_star = data->v_star;
    double** phi = data->phi;

    int i, j;

    for (i=1; i<m-1; i++){
        for (j=1; j<n-1; j++) {
            u_new[i][j] = u_star[i][j] - dt * (phi[i + 1][j] - phi[i][j]) / h;
            v_new[i][j] = v_star[i][j] - dt * (phi[i][j + 1] - phi[i][j]) / h;
            P[i][j] += phi[i][j];
        }
    }
    
    data->u_n_1 = u_n;
    data->u_n = u_new;
    data->v_n_1 = v_n;
    data->v_n = v_new;
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
