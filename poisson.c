//
// Created by Baptiste Hardy on 9/03/18.
//

#include "poisson.h"

//void old_poisson_solver(Data* data)
//{
//    int m = data->m;
//    int n = data->n;
//    double h = data->h;
//    double** u_star = data->u_star;
//    double** v_star = data->v_star;
//    double** phi = data->phi;
//
//    double **R = make2DDoubleArray(m,n);
//    double alpha = data->alpha;
//    double dt = data->dt;
//    double e;
//    int SORiter = 0;
//    double div;
//    double h2 = h*h;
//
//    do{
//        for (int i=1; i<m-1; i++){
//            for (int j=1; j<n-1; j++){
//                double phistar = .25*(-(h/dt)*(u_star[i][j]-u_star[i-1][j]+v_star[i][j]-v_star[i][j-1])
//                                 + phi[i+1][j] + phi[i-1][j] + phi[i][j+1] + phi[i][j-1]);
//
//                phi[i][j] = alpha*phistar + (1.-alpha)*phi[i][j];
//            }
//        }
//        //Adapt ghost points on Phi
//        for (int i=1; i<m-1; i++){
//            phi[i][0] = phi[i][1];
//            phi[i][n-1] = phi[i][n-2];
//        }
//
//        // LEFT AND RIGHT
//        for (int j=1; j<n-1; j++){
//            phi[0][j] = phi[1][j];
//            phi[m-1][j] = -phi[m-2][j];
//        }
//
//        double Sum = 0.;
//
//        for (int i=1; i<m-1; i++){
//            for (int j=1; j<n-1; j++){
//                  R[i][j] = -(u_star[i][j]-u_star[i-1][j]+v_star[i][j]-v_star[i][j-1])/(h*dt)
//                              +(phi[i+1][j]+phi[i-1][j]+phi[i][j+1]+phi[i][j-1]-4.*phi[i][j])/(h*h);
//
//                // CACLULER LA SOMME
//                Sum += R[i][j]*R[i][j];
//
//            }
//        }
//
//        e = (h*dt*d/Um)*pow(Sum/(L*d),0.5);
//        SORiter ++;
//
//        if(SORiter % 100 == 0)
//        {
//            printf("======= ======= iteration : %d done, Error : %f     \n", SORiter,e);
//        }
//
//
//    } while (e > data->SORtol && SORiter < data->SORitermax);
//    PetscPrintf(PETSC_COMM_WORLD, "ERROR  = %1.3e after %d/%d iterations \n",e, SORiter, data->SORitermax);

//}
