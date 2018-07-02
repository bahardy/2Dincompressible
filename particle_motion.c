//
// Created by Baptiste Hardy on 29/06/18.
//
#include "main.h"
#include "particle_motion.h"


void update_Up(Data* data, double* Up_k, double* Vp_k, double* Omega_p_k, int k)
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


    dudt[k] = (23.*F[k][2]-16.*F[k][1]+5.*F[k][0])/(12.*Sp[k]*(rho_r - 1.))
              + (23.*Fx_coll[k][2]-16.*Fx_coll[k][1]+5.*Fx_coll[k][0])/(12.*Sp[k]*(rho_p - rho_f)) - g;
    Up_k[k] = Up[k][0] + dt*dudt[k];
    dvdt[k] = (23.*G[k][2]-16.*G[k][1]+5.*G[k][0])/(12.*Sp[k]*(rho_r - 1.))
              + (23.*Fy_coll[k][2]-16.*Fy_coll[k][1]+5.*Fy_coll[k][0])/(12.*Sp[k]*(rho_p - rho_f)) ;
    Vp_k[k] = Vp[k][0] + dt*dvdt[k];
    domegadt[k] = (23.*M[k][2]-16.*M[k][1]+5.*M[k][0])/(12.*J[k]*(rho_r - 1.));
    Omega_p_k[k] = Omega_p[k][0] + dt*domegadt[k];
}

void update_Xp(Data* data, double* Xp_k, double* Yp_k, double* theta_k, int k)
{
    double* xg = data->xg;
    double* yg = data->yg;
    double* theta = data->theta;
    double** Up = data->Up;
    double** Vp = data->Vp;
    double** Omega_p = data->Omega_p;
    double dt = data->dt;

    Xp_k[k] = xg[k] + dt*(Up[k][1]+ Up[k][0])/2.;
    Yp_k[k] = yg[k] + dt*(Vp[k][1]+ Vp[k][0])/2.;
    theta_k[k] = theta[k] + dt*(Omega_p[k][1]+Omega_p[k][0])/2.;
    //PetscPrintf(PETSC_COMM_WORLD,"Position of the center of mass of particle %d: (x,y) = (%f,%f) \n", k+1, xg[k], yg[k]);
    //PetscPrintf(PETSC_COMM_WORLD,"Angle: theta  = %f \n", theta[k]);

}