//
// Created by Baptiste Hardy on 9/03/18.
//

#include "poisson.h"

void old_poisson_solver(Data* data, double** rho_new)
{
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
//    double div, drhodt;
//    double h2 = h*h;
//
//    do{
//        for (int i=1; i<m-1; i++){
//            for (int j=1; j<n-1; j++){
//                div = (.5*(rho_new[i][j]+rho_new[i+1][j])*u_star[i][j] - .5*(rho_new[i-1][j]+rho_new[i][j])*u_star[i-1][j])/h
//                           + (.5*(rho_new[i][j]+rho_new[i][j+1])*v_star[i][j] - .5*(rho_new[i][j-1]+rho_new[i][j])*v_star[i][j-1])/h;
//
//                double phistar = .25*h2*(-(div + drhodt)/dt
//                                      + (phi[i+1][j] + phi[i-1][j] + phi[i][j+1] + phi[i][j-1])/h2 );
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
//                div = (.5*(rho_new[i][j]+rho_new[i+1][j])*u_star[i][j] - .5*(rho_new[i-1][j]+rho_new[i][j])*u_star[i-1][j])/h
//                      + (.5*(rho_new[i][j]+rho_new[i][j+1])*v_star[i][j] - .5*(rho_new[i][j-1]+rho_new[i][j])*v_star[i][j-1])/h;
//                drhodt = (3*rho_new[i][j] - 4*rho_n[i][j] + rho_n_1[i][j])/(2.*dt);
//
//                R[i][j] = (phi[i+1][j]+phi[i-1][j]+phi[i][j+1]+phi[i][j-1]-4.*phi[i][j])/h2
//                        - (drhodt + div)/dt;
//                // CACLULER LA SOMME
//                Sum += R[i][j]*R[i][j];
//
//            }
//        }
//
//        e = (dt*data->d/(data->u_m * data->rho_f_ref))*sqrt(Sum*h2/(data->L*data->d));
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

}
