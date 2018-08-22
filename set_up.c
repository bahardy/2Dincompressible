//
// Created by Baptiste Hardy on 4/07/18.
//

#include "set_up.h"
#include "flow_solver.h"

void set_up(Data* data, int argc, char *argv[], int rank)
{
    /* DIMENSIONS */
    data->Dp = 1;
    data->d = 4.;
    data->H = 0.5*data->d;
    data->L = 6.;
    data->h = data->Dp/30;
    data->eps = 0;
#ifdef SMOOTHING
    data->eps = data->h;
#endif
    /* Collision paramaters */
    data->ep = 1e-3;
    data->Ep = 0.01*data->ep;
    data->ew = 1e-3; //data->h*data->h;//1e-6;
    data->Ew = 0.01*data->ew;//1e-8;

    /* PHYSICAL PARAMETERS */
    data->rho_f = 1.;
    data->rho_p = 100;
    data->rho_r = data->rho_p/data->rho_f;
    data->cp = 1000.;
    data->cf = 1000.;
    data->cr = data->cp/data->cf;

    // /* NON-DIMENSIONAL NUMBERS */
    data->Pr = 0.7;
    data->Le = 1; /* Lewis number, ratio between Sc and Prandtl */
    data->Sc = data->Le*data->Pr;
    data->Rep = 40.;
    data->Fr = sqrt(1e3);
    data->Da = 0.2;

    data->g =  1/pow((data->Fr),2);

#ifdef SEDIMENTATION
    data->Ga = 100*sqrt(data->rho_r -1);//39.15*sqrt(data->rho_r-1); //usually:100
    data->Rep = data->Ga;
    data->g = 1./(data->rho_r -1);
#endif

    /* FLOW */
    data->u_m = 1.;
    data->nu = data->u_m*data->Dp/data->Rep;

    /* ENERGY */
    data->alpha_f = data->nu/data->Pr;
    data->T0 = 1; //  temperature at the inlet
    data->Tp0 = 0; // initial particle temperature

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
    data->CA0 = 0.;
    data->CB0 = 0.;

    /* GRID */
    data->N = (int) (data->d/data->h);
    data->M = (int) (data->L/data->h);
    data->n = data->N + 2; /*for ghost points */
    data->m = data->M + 2; /*for ghost points */


    /* TIME INTEGRATION */
    data->CFL = 0.01; /*Courant-Freidrichs-Lewy condition on convective term */
    data->r = 0.2; /* Fourier condition on diffusive term */
    double dt_CFL = data->CFL*data->h/data->u_m;
    double dt_diff = data->r*data->h*data->h/data->nu;

#ifdef EXPLICIT
    data->ratio_dtau_dt = 1;
#endif
#ifndef EXPLICIT
    data->ratio_dtau_dt = 1e-4;
#endif
    data->dt = fmin(dt_CFL, dt_diff);
    data->dtau = data->ratio_dtau_dt*data->dt;

    if(rank == 0){
        printf("nu = %f\n", data->nu);
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
    data->t_move = 0; //data->Tf/
    data->t_coupling = 0;
    data->t_transfer = 0.;
    data->nKmax = 2;
    data->Kmax = 50; /* number of ramping steps */


    if(rank == 0){
        printf("Write every %d * dt \n", data->T_write);
        printf("Write %d times \n", data->N_write);
        printf("Final time : %f \n \n", data->Tf);
    }
}

void initialize_fields(Data* data)
{

    for (int k = 0; k < data->Np; k++) {
        data->dp[k] = data->Dp;
        data->rp[k] = .5 * data->Dp;

#ifdef DISK
        data->Sp[k] = M_PI * data->rp[k] * data->rp[k];
        //data->II[k]=(data->dp[k]*data->dp[k])/8.; /* IN 2-D !! */
        data->J[k] = (M_PI / 2.) * pow(data->rp[k], 4);
#endif
#ifdef ELLIPSE
        data->a = 2*Dp;
        data->b = Dp;
        data->Sp[k]=M_PI*data->a*data->b;
        data->II[k]=.25*(pow(data->a,2) + pow(data->b,2));
        data->J[k] =(M_PI/4.)*data->a*data->b*(pow(data->a,2) + pow(data->b,2));
#endif
    }

    /* Initialization of particles position */
    data->xg[0][0] = 2;
    data->xg[0][1] = data->xg[0][0];
    data->xg[0][2] = data->xg[0][1];

    data->yg[0][0] = data->H;
    data->yg[0][1] = data->yg[0][0];
    data->yg[0][2] = data->yg[0][1];

    data->theta[0][0] = 0;
    data->theta[0][1] = data->theta[0][0];
    data->theta[0][2] = data->theta[0][1];

//    data->xg[1][0] = 3;
//    data->xg[1][1] = data->xg[1][0];
//    data->xg[1][2] = data->xg[1][1];
//
//    data->yg[1][0] = data->H;
//    data->yg[1][1] = data->yg[1][0];
//    data->yg[1][2] = data->yg[1][1];
//
//    data->theta[1][0] = 0;
//    data->theta[1][1] = data->theta[1][0];
//    data->theta[1][2] = data->theta[1][1];
//
//    data->xg[2][0] = 7.5;
//    data->xg[2][1] = data->xg[2][0];
//    data->xg[2][2] = data->xg[2][1];
//
//    data->yg[2][0] = data->H;
//    data->yg[2][1] = data->yg[2][0];
//    data->yg[2][2] = data->yg[2][1];
//
//    data->theta[2][0] = 0;
//    data->theta[2][1] = data->theta[2][0];
//    data->theta[2][2] = data->theta[2][1];


    /*Initialization of particles velocities */
    for (int k = 0; k < data->Np; k++) {
        data->Up[k][2] = 0 * data->u_m;
        data->Up[k][1] = data->Up[k][2];
        data->Up[k][0] = data->Up[k][2];

        data->Vp[k][2] = 0 * data->u_m;
        data->Vp[k][1] = data->Vp[k][2];
        data->Vp[k][0] = data->Vp[k][2];

    }


    /* VELOCITY : horizontal flow Um  */
    for (int i = 0; i < data->m; i++) {
        for (int j = 0; j < data->n; j++) {
            data->u_n[i][j] = 1*data->u_m;
            data->u_n_1[i][j] = data->u_n[i][j];
            data->u_star[i][j] = data->u_n[i][j];
            data->T_n[i][j] = 0;
            data->T_n_1[i][j] = data->T_n[i][j];
            for (int s = 0; s < data->Ns; s++) {
                data->C_n[s][i][j] = 0;
                data->C_n_1[s][i][j] = data->C_n[s][i][j];
            }
            /* v_n is initially at zero */
        }
    }

    /*INLET BC*/
    for (int j = 0; j < data->n; j++) {
        if (j > (int) (0.5/data->h))
        {
            data->u_n[0][j] = 1*data->u_m;
            data->u_n_1[0][j] = data->u_n[0][j];
            data->u_star[0][j] = data->u_n[0][j];
        }
        else
        {
            data->u_n[0][j] = 0 * data->u_m;
            data->u_n_1[0][j] = data->u_n[0][j];
            data->u_star[0][j] = data->u_n[0][j];
        }
    }

    /*Initialization of particles temperatures */
    for (int k = 0; k < data->Np; k++) {
        data->Tp[k] = data->Tp0;
    }

    /* Inlet concentration */
    data->C0[0] = data->CA0;

}


