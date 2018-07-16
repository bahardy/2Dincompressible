//
// Created by Baptiste Hardy on 9/03/18.
//

#include "poisson.h"

void old_poisson_solver(Data* data)
{
    int m = data->m;
    int n = data->n;
    double d = data->d;
    double Um = data->u_m;
    double L = data->L;
    double h = data->h;
    double** u_star = data->u_star;
    double** v_star = data->v_star;
    double** phi = data->phi;

    double **R = make2DDoubleArray(m,n);
    double alpha = data->alpha_SOR;
    double dt = data->dt;
    double e;
    int SORiter = 0;
    double div;
    double h2 = h*h;

    do{
        for (int i=1; i<m-1; i++){
            for (int j=1; j<n-1; j++){
                double phistar = .25*(-(h/dt)*(u_star[i][j]-u_star[i-1][j]+v_star[i][j]-v_star[i][j-1])
                                 + phi[i+1][j] + phi[i-1][j] + phi[i][j+1] + phi[i][j-1]);

                phi[i][j] = alpha*phistar + (1.-alpha)*phi[i][j];
            }
        }
        //Adapt ghost points on Phi
        for (int i=1; i<m-1; i++){
            phi[i][0] = phi[i][1];
            phi[i][n-1] = phi[i][n-2];
        }

        // LEFT AND RIGHT
        for (int j=1; j<n-1; j++){
            phi[0][j] = phi[1][j];
            phi[m-1][j] = -phi[m-2][j];
        }

        double Sum = 0.;

        for (int i=1; i<m-1; i++){
            for (int j=1; j<n-1; j++){
                  R[i][j] = -(u_star[i][j]-u_star[i-1][j]+v_star[i][j]-v_star[i][j-1])/(h*dt)
                              +(phi[i+1][j]+phi[i-1][j]+phi[i][j+1]+phi[i][j-1]-4.*phi[i][j])/(h*h);

                // CACLULER LA SOMME
                Sum += R[i][j]*R[i][j];

            }
        }

        e = (h*dt*d/Um)*pow(Sum/(L*d),0.5);
        SORiter ++;

        if(SORiter % 100 == 0)
        {
            printf("======= ======= iteration : %d done, Error : %f     \n", SORiter,e);
        }


    } while (e > data->SORtol && SORiter < data->SORitermax);
    PetscPrintf(PETSC_COMM_WORLD, "ERROR  = %1.3e after %d/%d iterations \n",e, SORiter, data->SORitermax);
}

void old_poisson_solver_periodic(Data* data)
{
    int m = data->m;
    int n = data->n;
    double d = data->d;
    double Um = data->u_m;
    double L = data->L;
    double h = data->h;
    double** u_star = data->u_star;
    double** v_star = data->v_star;
    double** phi = data->phi;

    double **R = make2DDoubleArray(m,n);
    double alpha = data->alpha_SOR;
    double dt = data->dt;
    double e;
    int SORiter = 0;

    for (int i=0; i<m; i++) {
        for (int j = 1; j < n - 1; j++) {
            phi[i][j] = 0;
        }
    }
    do{
        for (int i=0; i<m; i++){
            for (int j=1; j<n-1; j++){
                double phistar = .25*(-(h/dt)*(u_star[i][j]-u_star[(i-1+m)%m][j]+v_star[i][j]-v_star[i][j-1])
                                      + phi[(i+1+m)%m][j] + phi[(i-1+m)%m][j] + phi[i][j+1] + phi[i][j-1]);

                phi[i][j] = alpha*phistar + (1.-alpha)*phi[i][j];
            }
        }
        //Adapt ghost points on Phi
        //along the walls
        for (int i=0; i<m; i++){
            phi[i][0] = phi[i][1];
            phi[i][n-1] = phi[i][n-2];
        }

        double Sum = 0.;

        for (int i=0; i<m; i++){
            for (int j=1; j<n-1; j++){
                R[i][j] = -(u_star[i][j]-u_star[(i-1+m)%m][j]+v_star[i][j]-v_star[i][j-1])/(h*dt)
                          +(phi[(i+1+m)%m][j]+phi[(i-1+m)%m][j]+phi[i][j+1]+phi[i][j-1]-4.*phi[i][j])/(h*h);
                // CACLULER LA SOMME
                Sum += R[i][j]*R[i][j];
            }
        }

        e = (h*dt*d/Um)*pow(Sum/(L*d),0.5);
        SORiter ++;

        if(SORiter % 100 == 0)
        {
            printf("======= ======= iteration : %d done, Error : %f     \n", SORiter,e);
        }


    } while (e > data->SORtol && SORiter < data->SORitermax);
    PetscPrintf(PETSC_COMM_WORLD, "ERROR  = %1.3e after %d/%d iterations \n",e, SORiter, data->SORitermax);
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
            phi[m-1][jj] = -phi[m-2][jj];
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

PetscErrorCode poisson_solver_periodic(Data* data, int myrank, int nbproc)
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

        MatSetValue(A, r, (r-N + M*N)%(M*N), -1., INSERT_VALUES);
        if(jj>0){
            MatSetValue(A, r, r-1, -1., INSERT_VALUES);
        }
        MatSetValue(A, r, r, 4., INSERT_VALUES);
        if(jj<N-1){
            MatSetValue(A, r, r+1, -1., INSERT_VALUES);
        }
        MatSetValue(A, r, (r+N + M*N)%(M*N), -1., INSERT_VALUES);
        if(jj == 0 || jj == N-1){
            MatSetValue(A, r, r, 3., INSERT_VALUES);
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
        i = ii; // periodic case
        j = jj+1;
        div_u_star = (u_star[i][j]-u_star[(i-1+m)%m][j]+v_star[i][j]-v_star[i][j-1])/h;

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
            i = r/N;
            jj = r%N; j = jj+1;
            phi[i][j] = array[r];
        }
        // update ghost points on phi
        for(i=0; i<m; i++){
            /* cancel gradients : dp/dn=0 --> dphi/dn = 0*/
            phi[i][0] = phi[i][1];
            phi[i][n-1] = phi[i][n-2];
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

        for(r = 0; r < M*N; r++){
            i = r/N; j = r%N;
            jj = j+1;
            phi[i][jj] = array[r];
        }
        // update ghost points on phi
        for(i=0; i<m; i++){
            /* cancel gradients : dp/dn=0 --> dphi/dn = 0*/
            phi[i][0] = phi[i][1];
            phi[i][n-1] = phi[i][n-2];
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

void poisson_residual(Data* data)
{
    double h = data->h;
    double** u_star = data->u_star;
    double** v_star = data->v_star;
    double** phi = data->phi;
    double dt = data->dt;
    double d = data->d;
    double L = data->L;
    double u_m = data->u_m;

    int n = data->n;
    int m = data->m;
    double** res = make2DDoubleArray(m,n);
    double e = 0;
    double sum = 0.;

    for (int i=1; i<m-1; i++){
        for (int j=1; j<n-1; j++){
            res[i][j] = (phi[i+1][j]+phi[i-1][j]+phi[i][j+1]+phi[i][j-1]-4.*phi[i][j])/(h*h)
                      -(u_star[i][j]-u_star[i-1][j]+v_star[i][j]-v_star[i][j-1])/(h*dt);

            sum += res[i][j]*res[i][j];

        }
    }

    e = (h*dt*d/u_m)*sqrt(sum/(L*d));

    PetscPrintf(PETSC_COMM_WORLD, "Relative residual e = %1e6 \n", e);
    free2Darray(res, m);

}

void poisson_residual_periodic(Data* data)
{
    double h = data->h;
    double** u_star = data->u_star;
    double** v_star = data->v_star;
    double** phi = data->phi;
    double dt = data->dt;
    double d = data->d;
    double L = data->L;
    double u_m = data->u_m;

    int n = data->n;
    int m = data->m;
    double** res = make2DDoubleArray(m,n);
    double e = 0;
    double sum = 0.;

    for (int i=0; i<m; i++){
        for (int j=1; j<n-1; j++){
            res[i][j] = (phi[(i+1+m)%m][j]+phi[(i-1+m)%m][j]+phi[i][j+1]+phi[i][j-1]-4.*phi[i][j])/(h*h)
                        -(1./dt)*(u_star[i][j]-u_star[(i-1+m)%m][j]+v_star[i][j]-v_star[i][j-1])/(h);

            sum += res[i][j]*res[i][j];

        }
    }

    //e = (dt*d/u_m)*sqrt(sum*h*h/(L*d));
    e = sqrt(sum*h*h/(L*d));
    PetscPrintf(PETSC_COMM_WORLD, "Relative residual e = %1e6 \n", e);
    free2Darray(res, m);

}