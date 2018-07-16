//
// Created by Baptiste Hardy on 6/07/18.
//

#include "set_up_periodic.h"

void set_up_periodic(Data *data, int argc, char **argv, int rank)
{
    /* DIMENSIONS */
    data->Dp = 1;
    data->d = 4.;
    data->H = 0.5*data->d;
    data->L = 8.;
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
    data->rho_p = 1.5;
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
    data->CFL = 0.05;  /*Courant-Freidrichs-Lewy condition on convective term */
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
    //data->dt = 0.005;
    data->dtau = data->ratio_dtau_dt*data->dt;

    data->SORitermax = 100000;
    data->alpha_SOR = 1.98;
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

void initialize_fields_periodic(Data* data)
{
    /* Particles position */
    data->xg[0][0] = 3;
    data->xg[0][1] = data->xg[0][0];
    data->xg[0][2] = data->xg[0][1];

    //data->xg[1] = 4.5;
    data->yg[0][0] = data->H;
    data->yg[0][1] = data->yg[0][0];
    data->yg[0][2] = data->yg[0][1];

    data->theta[0][0] = 0;
    data->theta[0][1] = data->theta[0][0];
    data->theta[0][2] = data->theta[0][1];

    data->dp[0] = data->Dp;
    //data->dp[1] = data->Dp;
    data->rp[0] = .5*data->Dp;
    //data->rp[1] = .5*data->Dp;

    data->Up[0][1] = 0.5*data->u_m;
    data->Up[0][0] = data->Up[0][1];
    data->Up[0][2] = data->Up[0][1];

    data->Vp[0][1] = data->u_m;
    data->Vp[0][0] = data->Vp[0][1];
    data->Vp[0][2] = data->Vp[0][1];


#ifdef TEMP
    /*Initialization of particles temperatures */
    for(int k=0; k<data->Np; k++){
        data->Tp[k] = data->Tp0;
    }

    /* We feed reactants at the inlet */
    data->C0[0] = data->CA0;
    data->C0[1] = data->CB0;
#endif

    for(int k=0; k<data->Np; k++){
#ifdef DISK
        data->Sp[k]=M_PI*data->rp[k]*data->rp[k];
        data->J[k] =(M_PI/2.)*pow(data->rp[k],4);
#endif
#ifdef ELLIPSE
        data->a = 2*Dp;
        data->b = Dp;
        data->Sp[k]=M_PI*data->a*data->b;
        data->J[k] =(M_PI/4.)*data->a*data->b*(pow(data->a,2) + pow(data->b,2));
#endif

    }
}