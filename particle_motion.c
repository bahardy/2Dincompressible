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

    if (data->iter <= 2)
    {
        xg[k][2] = xg[k][1] + dt*Up[k][2];
        yg[k][2] = yg[k][1] + dt*Vp[k][2];
        theta[k][2] = theta[k][1] + dt*Omega_p[k][2];
    }
    else {

#ifdef ITERATIVE
        xg[k][2] = xg[k][1] + dt*Up[k][2];
        yg[k][2] = yg[k][1] + dt*Vp[k][2];
        theta[k][2] = theta[k][1] + dt+Omega_p[k][2];
#endif

#ifdef LF
        xg[k][2] = xg[k][0] + 2*dt*Up[k][1];
        yg[k][2] = yg[k][0] + 2*dt*Vp[k][1];
        theta[k][2] = theta[k][0] + 2*dt+Omega_p[k][1];

    //    xg[k][2] = xg[k][1] + dt*Up[k][2];
    //    yg[k][2] = yg[k][1] + dt*Vp[k][2];
    //    theta[k][2] = theta[k][1] + dt*Omega_p[k][2];
#endif

#ifdef AB3
        xg[k][2] = xg[k][1] + dt * (23. * Up[k][2] - 16. * Up[k][1] + 5. * Up[k][0]) / 12.;
        yg[k][2] = yg[k][1] + dt * (23. * Vp[k][2] - 16. * Vp[k][1] + 5. * Vp[k][0]) / 12.;
        theta[k][2] = theta[k][1] + dt * (23. * Omega_p[k][2] - 16. * Omega_p[k][1] + 5. * Omega_p[k][0]) / 12.;

#endif

#ifdef EE
        xg[k][2] = xg[k][1] + dt*Up[k][2];
        yg[k][2] = yg[k][1] + dt*Vp[k][2];
        theta[k][2] = theta[k][1] + dt*Omega_p[k][2];
#endif

    }
}

void update_Up(Data* data, int k)
{
    double rho_r = data->rho_r;
    double rho_f = data->rho_f;
    double rho_p = data->rho_s;
    double* Sp = data->Sp;
    double* J = data->J; //m^4
    double** F = data->F;
    double** G = data->G;
    double** M = data->Mz;
#ifdef GRAVITY
    double g = data->g;
#else
    double g = 0;
#endif

    double** Up = data->Up;
    double** Vp = data->Vp;
    double** Omega_p = data->Omega_p;
    double** Fx_coll = data->Fx_coll;
    double** Fy_coll = data->Fy_coll;
    double dudt, dvdt, domegadt;
    double** ax = data->dudt;
    double** ay = data->dvdt;
    double** atheta = data->domegadt;

    double dt = data->dt;

#ifndef ITERATIVE
#ifndef LF    
    data->Up[k][0] = data->Up[k][1];
    data->Up[k][1] = data->Up[k][2];

    data->Vp[k][0] = data->Vp[k][1];
    data->Vp[k][1] = data->Vp[k][2];

    data->Omega_p[k][0] = data->Omega_p[k][1];
    data->Omega_p[k][1] = data->Omega_p[k][2];
#endif
#endif

#ifdef ITERATIVE
    dudt = F[k][2]/(Sp[k]*(rho_r-1)) + Fx_coll[k][2]/(Sp[k]*(rho_s - rho_f)) -g;
    dvdt = G[k][2]/(Sp[k]*(rho_r-1)) + Fy_coll[k][2]/(Sp[k]*(rho_s - rho_f));
    domegadt = M[k][2]/(J[k]*(rho_r-1));

    Up[k][2] = Up[k][1] + dt*dudt;
    Vp[k][2] = Vp[k][1] + dt*dvdt;
    Omega_p[k][2] = Omega_p[k][1] + dt*domegadt;
#else
    if (data->iter <= 2)
    {
        dudt = F[k][2]/(Sp[k]*(rho_r-1)) + Fx_coll[k][2]/(Sp[k]*(rho_p - rho_f)) - g;
        dvdt = G[k][2]/(Sp[k]*(rho_r-1)) + Fy_coll[k][2]/(Sp[k]*(rho_p - rho_f));
        domegadt = M[k][2]/(J[k]*(rho_r-1));

        Up[k][2] = Up[k][1] + dt*dudt;
        Vp[k][2] = Vp[k][1] + dt*dvdt;
        Omega_p[k][2] = Omega_p[k][1] + dt*domegadt;
    }
    else
    {
#ifdef LF
        dudt = F[k][2]/(Sp[k]*(rho_r-1)) + Fx_coll[k][2]/(Sp[k]*(rho_s - rho_f));
        dvdt = G[k][2]/(Sp[k]*(rho_r-1)) + Fy_coll[k][2]/(Sp[k]*(rho_s - rho_f));
        domegadt = M[k][2]/(J[k]*(rho_r-1));

        Up[k][2] = Up[k][0] + 2*dt*dudt;
        Vp[k][2] = Vp[k][0] + 2*dt*dvdt;
        Omega_p[k][2] = Omega_p[k][0] + 2*dt*domegadt;
#endif

#ifdef AB3
        dudt = (23.*F[k][2]-16.*F[k][1]+5.*F[k][0])/(12.*Sp[k]*(rho_r - 1.)) +
                (23.*Fx_coll[k][2]-16.*Fx_coll[k][1]+5.*Fx_coll[k][0])/(12.*Sp[k]*(rho_p - rho_f)) - g;
        dvdt = (23.*G[k][2]-16.*G[k][1]+5.*G[k][0])/(12.*Sp[k]*(rho_r - 1.)) +
                (23.*Fy_coll[k][2]-16.*Fy_coll[k][1]+5.*Fy_coll[k][0])/(12.*Sp[k]*(rho_p - rho_f));
        domegadt = (23.*M[k][2]-16.*M[k][1]+5.*M[k][0])/(12.*J[k]*(rho_r - 1.));

//        dudt = (1./rho_r)*ax[k][2]
//               + (23.*F[k][2]-16.*F[k][1]+5.*F[k][0])/(12.*Sp[k]*rho_r)
//               + (23.*Fx_coll[k][2]-16.*Fx_coll[k][1]+5.*Fx_coll[k][0])/(12.*Sp[k]*rho_s)
//               - (1 - 1/rho_r)*g;
//
//        dvdt = (1./rho_r)*ay[k][2]
//               + (23.*G[k][2]-16.*G[k][1]+5.*G[k][0])/(12.*Sp[k]*rho_r)
//               + (23.*Fy_coll[k][2]-16.*Fy_coll[k][1]+5.*Fy_coll[k][0])/(12.*Sp[k]*rho_s);
//
//        domegadt = (1./rho_r)*atheta[k][2]
//                   + (23.*M[k][2]-16.*M[k][1]+5.*M[k][0])/(12.*J[k]*rho_r);

        Up[k][2] = Up[k][1] + dt*dudt;
        Vp[k][2] = Vp[k][1] + dt*dvdt;
        Omega_p[k][2] = Omega_p[k][1] + dt*domegadt;
#endif

#ifdef EE
//        dudt = (1./rho_r)*ax[k][2] + F[k][2]/(Sp[k]*rho_r) + Fx_coll[k][2]/(Sp[k]*rho_s)
//               -(1 - 1./rho_r)*g;
//
//        dvdt = (1./rho_r)*ay[k][2] + G[k][2]/(Sp[k]*rho_r) + Fy_coll[k][2]/(Sp[k]*rho_s);
//
//        domegadt = (1./rho_r)*aomega[k][2] + M[k][2]/(J[k]*rho_r);

	dudt = F[k][2]/(Sp[k]*(rho_r-1)) + Fx_coll[k][2]/(Sp[k]*(rho_s - rho_f));
	dvdt = G[k][2]/(Sp[k]*(rho_r-1)) + Fy_coll[k][2]/(Sp[k]*(rho_s - rho_f));
	domegadt = M[k][2]/(J[k]*(rho_r-1));

        Up[k][2] = Up[k][1] + dt*dudt;
        Vp[k][2] = Vp[k][1] + dt*dvdt;
        Omega_p[k][2] = Omega_p[k][1] + dt*domegadt;
#endif

    }
#endif

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
    double rho_p = data->rho_s;
    double rho_f = data->rho_f;
    double cr = data->cr;
    double cp = data->cp_s;
    double cf = data->cp_f;
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

    //dCdt[k][1] = (23.*(PP[k][0][2]+PP[k][1][2])-16.*(PP[k][0][1]+PP[k][1][1])+5.*(PP[k][0][0]+PP[k][1][0]))/12.;
    Cp[k][0] = 0.;
    //Cp[k][1] += dt*dCdt[k][1];
}

void compute_Qr(double** Qr, double rate, double dH, int k)
{

    Qr[k][0] = Qr[k][1];
    Qr[k][1] = Qr[k][2];
    Qr[k][2] = rate*(-dH);
}
