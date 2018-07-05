// Created by Baptiste Hardy on 29/06/18.
//
#include "main.h"
#include "particle_motion.h"



void update_Xp(Data* data, int k)
{
    double** xg = data->xg;
    double** yg = data->yg;
    double** theta = data->theta;
    double** Up = data->Up;
    double** Vp = data->Vp;
    double** Omega_p = data->Omega_p;
    double dt = data->dt;

    xg[k][1] = xg[k][0] + dt*(Up[k][2]+ Up[k][1])/2.;
    yg[k][1] = yg[k][0] + dt*(Vp[k][2]+ Vp[k][1])/2.;
    theta[k][1] = theta[k][0] + dt*(Omega_p[k][2]+Omega_p[k][1])/2.;

    /*xg[k] += dt*(23.*Up[k][2]-16.*Up[k][1]+5.*Up[k][0])/12.;
    yg[k] += dt*(23.*Vp[k][2]-16.*Vp[k][1]+5.*Vp[k][0])/12.;
    theta[k] += dt*(23.*Omega_p[k][2]-16.*Omega_p[k][1]+5.*Omega_p[k][0])/12.;*/

    //PetscPrintf(PETSC_COMM_WORLD,"Position of the center of mass of particle %d: (x,y) = (%f,%f) \n", k+1, xg[k], yg[k]);
    //PetscPrintf(PETSC_COMM_WORLD,"Angle: theta  = %f \n", theta[k]);

}

void update_Up(Data* data, int k)
{
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
    double dudt, dvdt, domegadt;
    double dt = data->dt;

    dudt = F[k][1]/(Sp[k]*(rho_r-1)) - g + Fx_coll[k][1]/(Sp[k]*(rho_p - rho_f));
    dvdt = G[k][1]/(Sp[k]*(rho_r-1)) + Fy_coll[k][1]/(Sp[k]*(rho_p - rho_f));
    domegadt = M[k][1]/(J[k]*(rho_r-1));

#ifdef ITERATIVE
    Up[k][2] = Up[k][1] + dt*dudt;
    Vp[k][2] = Vp[k][1] + dt*dvdt;
    Omega_p[k][2] = Omega_p[k][1] + dt*domegadt;
#endif

#ifndef ITERATIVE
    if (data->iter == 1) {
        Up[k][2] = Up[k][1] + dt*dudt;
        Vp[k][2] = Vp[k][1] + dt*dvdt;
        Omega_p[k][2] = Omega_p[k][1] + dt*domegadt;
    }
    else
    {
        Up[k][2] = Up[k][0] + 2*dt*dudt;
        Vp[k][2] = Vp[k][0] + 2*dt*dvdt;
        Omega_p[k][2] = Omega_p[k][0] + 2*dt*domegadt;
    }
#endif

    /* dudt[k] = (23.*F[k][2]-16.*F[k][1]+5.*F[k][0])/(12.*Sp[k]*(rho_r - 1.)) + (23.*Fx_coll[k][2]-16.*Fx_coll[k][1]+5.*Fx_coll[k][0])/(12.*Sp[k]*(rho_p - rho_f)) - g;
     Up[k][3] = Up[k][2] + dt*dudt[k];
     dvdt[k] = (23.*G[k][2]-16.*G[k][1]+5.*G[k][0])/(12.*Sp[k]*(rho_r - 1.)) + (23.*Fy_coll[k][2]-16.*Fy_coll[k][1]+5.*Fy_coll[k][0])/(12.*Sp[k]*(rho_p - rho_f)) ;
     Vp[k][3] = Vp[k][2] + dt*dvdt[k];
     domegadt[k] = (23.*M[k][2]-16.*M[k][1]+5.*M[k][0])/(12.*J[k]*(rho_r - 1.));
     Omega_p[k][3] = Omega_p[k][2] + dt*domegadt[k];

     Up[k][0] = Up[k][1]; Up[k][1] = Up[k][2]; Up[k][2] = Up[k][3];
     Vp[k][0] = Vp[k][1]; Vp[k][1] = Vp[k][2]; Vp[k][2] = Vp[k][3];
     Omega_p[k][0] = Omega_p[k][1]; Omega_p[k][1] = Omega_p[k][2]; Omega_p[k][2] = Omega_p[k][3];*/
}

void relax_Up(Data* data, double relax, double* Up_old, double* Vp_old, double* Omega_p_old, int k)
{
    data->Up[k][2] = relax*data->Up[k][2] + (1-relax)*Up_old[k];
    data->Vp[k][2] = relax*data->Vp[k][2] + (1-relax)*Vp_old[k];
    data->Omega_p[k][2] = relax*data->Omega_p[k][2] + (1-relax)*Omega_p_old[k];
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

void compute_Qr(double** Qr, double rate, double dH, int k)
{

    Qr[k][0] = Qr[k][1];
    Qr[k][1] = Qr[k][2];
    Qr[k][2] = rate*(-dH);
}
