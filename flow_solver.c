//
// Created by Baptiste Hardy on 5/07/18.
//

#include "main.h"
#include "flow_solver.h"

void get_ghosts(Data* data, double T0, double* C0)
{
    double** u_n = data->u_n;
    double** v_n = data->v_n;
    double** T_n = data->T_n;
    double*** C = data->C_n;

    int m = data->m;
    int n = data->n;
    int Ns = data->Ns;

    /*Ghost points to impose BC's */
    /* Along y = -H and y = H */
    for (int i=0; i<m; i++){
        /* On u_n */
#ifdef SLIP
        /* Bottom and top Wall : slip : du/dn= 0 */
        u_n[i][0] = u_n[i][1];
        u_n[i][n-1] = u_n[i][n-2];
#endif
#ifndef SLIP
        /* No-slip (channel) */
        u_n[i][0] = -0.2*(u_n[i][3] - 5.*u_n[i][2] + 15.*u_n[i][1]);
        u_n[i][n-1] = -0.2*(u_n[i][n-4] - 5.*u_n[i][n-3] + 15.*u_n[i][n-2]);
#endif

#ifdef TEMP
        /* Walls : adiabatic: dTdn = 0, no mass flux */
        T_n[i][0] = T_n[i][1];
        T_n[i][n-1] = T_n[i][n-2];

        for (int s=0; s<Ns; s++){
            C[s][i][0] = C[s][i][1];
            C[s][i][n-1] = C[s][i][n-2];
        }
#endif
    }

    /* Along x = 0 and x = L */
    for (int j = 0; j<n; j++){
        /* Inflow : horizontal flow --> v_n = 0 */
        v_n[0][j] = -0.2*(v_n[3][j] - 5.*v_n[2][j] + 15.*v_n[1][j]);
        /* Natural outflow : dV/dx =0 */
        v_n[m-1][j] = v_n[m-2][j];

#ifdef TEMP
        /* On T_n and C */
        /* Inflow : T_n uniform  */
        T_n[0][j] = -0.2*(T_n[3][j]-5.*T_n[2][j]+15.*T_n[1][j]-16.*T0);
        T_n[m-1][j] = (7.*T_n[m-2][j]-5.*T_n[m-3][j]+T_n[m-4][j])/3.;

//        /* Inflow : CA = CA0; CB = CB0 */
//        C[0][0][j] = -0.2*(C[0][3][j]-5.*C[0][2][j]+15.*C[0][1][j]-16.*C0[0]);
//        C[1][0][j] = -0.2*(C[1][3][j]-5.*C[1][2][j]+15.*C[1][1][j]-16.*C0[1]);

        /*Outflow : We cancel axial dispersion d2T/dx2 = 0; d2C/dx2 = 0; */
        for(int s=0; s<Ns; s++){
            C[s][0][j] = -0.2*(C[s][3][j]-5.*C[s][2][j]+15.*C[s][1][j]-16.*C0[s]);
            C[s][m-1][j] = (7.*C[s][m-2][j]-5.*C[s][m-3][j]+C[s][m-4][j])/3.;
        }
#endif

    }
}

void get_masks(Data* data)
{
    double*** chi_S = data->chi_S;
    double*** chi_U = data->chi_U;
    double*** chi_V = data->chi_V;

    double** I_S = data->I_S;
    double** I_U = data->I_U;
    double** I_V = data->I_V;
    double*** Ip_S = data->Ip_S;
    double*** Ip_U = data->Ip_U;
    double*** Ip_V = data->Ip_V;

    double** coloring = data->coloring;
    double** xg = data->xg;
    double** yg = data->yg;
    double** theta = data->theta;
    double* rp = data->rp;
    double d;
    int m = data->m;
    int n = data->n;
    double h = data->h;
    int Np = data->Np;

    double xU, xV, yU, yV, yS, xS;
    double xG, yG;
    double xloc, yloc, delta;

    for(int i=0; i<m; i++){
        xU = i*h;
        xV = (i-0.5)*h;
        for(int j=0; j<n; j++){
            yU = (j-0.5)*h;
            yV = j*h;
            yS = yU;
            xS = xV;

            I_S[i][j] = 0; /*Reset the masks */
            I_U[i][j] = 0;
            I_V[i][j] = 0;
            coloring[i][j] = 0;

            /*Go over all the particles */
            for(int k=0; k<Np; k++){
                xG = xg[k][1];
                yG = yg[k][1];

#ifdef ELLIPSE
                double x;
		double y;
		/*ELLIPSE*/
		x = xU-xG;
		y = yU-yG;
                double EU = pow(b,2)*(pow(y*cos(theta[k]),2)+pow(x*sin(theta[k]),2)-y*x*sin(2*theta[k]))
                + pow(a,2)*(pow(x*cos(theta[k]),2)+pow(y*sin(theta[k]),2)+y*x*sin(2*theta[k]))
                - pow(a,2)*pow(b,2);

		x = xV-xG;
		y = yV-yG;
                double EV = pow(b,2)*(pow(y*cos(theta[k]),2)+pow(x*sin(theta[k]),2)-y*x*sin(2*theta[k]))
                + pow(a,2)*(pow(x*cos(theta[k]),2)+pow(y*sin(theta[k]),2)+y*x*sin(2*theta[k]))
                - pow(a,2)*pow(b,2);

		x = xS-xG;
		y = yS-yG;
                double ES = pow(b,2)*(pow(y*cos(theta[k]),2)+pow(x*sin(theta[k]),2)-y*x*sin(2*theta[k]))
                + pow(a,2)*(pow(x*cos(theta[k]),2)+pow(y*sin(theta[k]),2)+y*x*sin(2*theta[k]))
                - pow(a,2)*pow(b,2);

                Ip_S[k][i][j]=(ES <= 0);
                Ip_U[k][i][j]=(EU <= 0);
                Ip_V[k][i][j]=(EV <= 0);

		coloring[i][j] +=Ip_S[k][i][j];
#endif
#ifdef DISK
                chi_S[k][i][j]=((xS-xG)*(xS-xG)+(yS-yG)*(yS-yG)<= (rp[k]+data->eps)*(rp[k]+data->eps));
                chi_U[k][i][j]=((xU-xG)*(xU-xG)+(yU-yG)*(yU-yG)<= (rp[k]+data->eps)*(rp[k]+data->eps));
                chi_V[k][i][j]=((xV-xG)*(xV-xG)+(yV-yG)*(yV-yG)<= (rp[k]+data->eps)*(rp[k]+data->eps));

#ifndef SMOOTHING
                Ip_S[k][i][j]=((xS-xG)*(xS-xG)+(yS-yG)*(yS-yG)<= rp[k]*rp[k]);
                Ip_U[k][i][j]=((xU-xG)*(xU-xG)+(yU-yG)*(yU-yG)<= rp[k]*rp[k]);
                Ip_V[k][i][j]=((xV-xG)*(xV-xG)+(yV-yG)*(yV-yG)<= rp[k]*rp[k]);
#endif

#ifdef SMOOTHING

                //Smoothing S
                d = rp[k] - sqrt((xS-xG)*(xS-xG)+(yS-yG)*(yS-yG));
                if( d < - data->eps)
                    Ip_S[k][i][j] = 0;
                else if( fabs(d) <= data->eps)
                    Ip_S[k][i][j] = .5*(1 + d/data->eps + (1./M_PI)*sin( M_PI* d/data->eps) );
                else if( d > data->eps)
                    Ip_S[k][i][j] = 1;

                //Smoothing U
                d = rp[k] - sqrt((xU-xG)*(xU-xG)+(yU-yG)*(yU-yG));
                if( d < - data->eps)
                    Ip_U[k][i][j] = 0;
                else if( fabs(d) <= data->eps)
                    Ip_U[k][i][j] = .5*(1 + d/data->eps + (1./M_PI)*sin( M_PI* d/data->eps) );
                else if( d > data->eps)
                    Ip_U[k][i][j] = 1;

                //Smoothing V
                d = rp[k] - sqrt((xV-xG)*(xV-xG)+(yV-yG)*(yV-yG));
                if( d < - data->eps)
                    Ip_V[k][i][j] = 0;
                else if( fabs(d) <= data->eps)
                    Ip_V[k][i][j] = .5*(1 + d/data->eps + (1./M_PI)*sin( M_PI* d/data->eps) );
                else if( d > data->eps)
                    Ip_V[k][i][j] = 1;

#endif

                xloc = xS-xG;
                yloc = yS-yG;
                delta = atan2(yloc, xloc);
                coloring[i][j] += Ip_S[k][i][j];

                if((int) floor((delta-theta[k][1])/(M_PI/2.)) % 2 == 0 ){
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

void get_Cs(Data* data)
{
    double*** Cs = data-> Cs;
    double** Cp = data->Cp;
    double*** chi_S = data->chi_S;
    int m = data->m;
    int n = data->n;
    int Np = data->Np;

    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            Cs[1][i][j] = 0.;
            for(int k=0; k<Np; k++){
                Cs[1][i][j] += chi_S[k][i][j]*Cp[k][1];
            }
        }
    }

}

void get_Ts(Data* data)
{
    double** Ts = data-> Ts;
    double* Tp = data->Tp;
    double*** chi_S = data->chi_S;
    int m = data->m;
    int n = data->n;
    int Np = data->Np;

    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            Ts[i][j] = 0.;
            for(int k = 0; k<Np; k++){
                Ts[i][j]+= chi_S[k][i][j]*Tp[k];

            }
        }
    }
}

void get_Us_Vs(Data* data)
{

    double** u_s = data-> u_s;
    double** v_s = data-> v_s;
    double** Up = data->Up;
    double** Vp = data->Vp;
    double** Omega_p = data->Omega_p;
    double** xg = data->xg;
    double** yg = data->yg;

    double*** chi_U = data->chi_U;
    double*** chi_V = data->chi_V;

    int m = data->m;
    int n = data->n;
    int Np = data->Np;
    double h = data->h;

    double xV, yU;
    for(int i=0; i<m; i++){
        xV = (i-0.5)*h;
        for(int j=0; j<n; j++){
            yU = (j-0.5)*h;
            u_s[i][j] = 0.;
            v_s[i][j] = 0.;
            for (int k = 0; k<Np; k++){
                u_s[i][j]+= chi_U[k][i][j]*(Up[k][3] - Omega_p[k][3]*(yU-yg[k][1]));
                v_s[i][j]+= chi_V[k][i][j]*(Vp[k][3] + Omega_p[k][3]*(xV-xg[k][1]));
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
    double nu = data->nu;

    double** I_U = data->I_U;
    double** I_V = data->I_V;

    double** u_n = data->u_n;
    double** u_star = data->u_star;
    double** u_s = data->u_s;

    double** v_n = data->v_n;
    double** v_star = data->v_star;
    double** v_s = data->v_s;

    double** P = data->P;

    double Um = data->u_m;

    double H_U, H_U_old;
    double H_V, H_V_old;
    double lapU, lapV;
    double dpdx, dpdy;

    double uL, uR, vT, vB;
    double dudxR, dudxL, dudyT, dudyB, dvdxR, dvdxL, dvdyT, dvdyB;

    /* u_star ADAMS-BASHFORTH 2 */
    for (i=1; i<m-2; i++){
        for (j=1; j<n-1; j++){

            uR = .5*(u_n[i][j] + u_n[i+1][j]);
            uL = .5*(u_n[i][j] + u_n[i-1][j]);
            dudxR = (u_n[i+1][j]- u_n[i][j])/h;
            dudxL = (u_n[i][j]- u_n[i-1][j])/h;

            vT = .5*(v_n[i][j] + v_n[i+1][j]);
            vB = .5*(v_n[i][j-1] + v_n[i+1][j-1]);
            dudyT = (u_n[i][j+1] - u_n[i][j])/h;
            dudyB = (u_n[i][j] - u_n[i][j-1])/h;

            H_U = .5*(uR*dudxR + uL*dudxL) + .5*(vT*dudyT + vB*dudyB);

            if (data->iter == 1){
                H_U_old = H_U;
            }
            else
            {
                H_U_old = data->H_u_n_1[i][j];
            }

            // LAPLACIAN
            lapU = (u_n[i+1][j]+u_n[i-1][j]+u_n[i][j+1]+u_n[i][j-1]-4.*u_n[i][j])/(h*h);

            // PRESSURE term
            dpdx = (P[i+1][j]-P[i][j])/h;

#ifdef EXPLICIT
            // EXPLICIT VERSION
            u_star[i][j] = u_n[i][j] + dt*(-1.5*H_U + 0.5*H_U_old - dpdx + nu*lapU) - ramp*I_U[i][j]*(dt/dtau)*(u_n[i][j] - u_s[i][j]);
#else
            // IMPLICIT VERSION
            u_star[i][j] = (u_n[i][j] + dt*(-1.5*H_U + 0.5*H_U_old - dpdx + nu*lapU) + (dt/dtau)*ramp*I_U[i][j]*u_s[i][j])/(1.+ramp*I_U[i][j]*dt/dtau);
#endif
            data->H_u_n_1[i][j] = H_U;
        }
    }

    /*Outflow condition */
    for (j=1; j<n-1; j++){
        u_star[m-2][j] = u_n[m-2][j] - dt*Um*(u_n[m-2][j]-u_n[m-3][j])/h;
    }
    /*u_star[0][j] (inflow) is fixed once for all at the beginning */

    /* v_star  ADAMS-BASHFORTH 2 */
    for (i=1; i<m-1; i++){
        for (j=1; j<n-2; j++){

            uR = .5*(u_n[i][j] + u_n[i][j+1]);
            uL = .5*(u_n[i-1][j] + u_n[i-1][j+1]);
            dvdxR = (v_n[i+1][j]- v_n[i][j])/h;
            dvdxL = (v_n[i][j]- v_n[i-1][j])/h;

            vT = .5*(v_n[i][j] + v_n[i][j+1]);
            vB = .5*(v_n[i][j] + v_n[i][j-1]);
            dvdyT = (v_n[i][j+1] - v_n[i][j])/h;
            dvdyB = (v_n[i][j] - v_n[i][j-1])/h;

            H_V = .5*(uR*dvdxR + uL*dvdxL) + .5*(vT*dvdyT + vB*dvdyB);

            if (data->iter == 1){
                H_V_old = H_V;
            }
            else
            {
                H_V_old = data->H_v_n_1[i][j];
            }

            // LAPLACIAN
            lapV = (v_n[i+1][j]+v_n[i-1][j]+v_n[i][j+1]+v_n[i][j-1]-4.*v_n[i][j])/(h*h);

            // PRESSURE TERM
            dpdy = (P[i][j+1]-P[i][j])/h;

#ifdef EXPLICIT
            // EXPLICIT VERSION
            v_star[i][j] = v_n[i][j] + dt*(-1.5*H_V + 0.5*H_V_old - dpdy + nu*lapV) - ramp*I_V[i][j]*(dt/dtau)*(v_n[i][j] - v_s[i][j]);
#else
            // IMPLICIT VERSION
            v_star[i][j] = (v_n[i][j] + dt*(-1.5*H_V + 0.5*H_V_old - dpdy + nu*lapV) + (dt/dtau)*ramp*I_V[i][j]*v_s[i][j])/(1.+ramp*I_V[i][j]*dt/dtau);
#endif
            /* the value of v_star on the boundaries (j=0, j=n-2) is set to zero at allocation */

            data->H_v_n_1[i][j] = H_V;
        }
    }
}

void update_flow(Data* data)
{

    int m = data->m;
    int n = data->n;
    double dt = data->dt;
    double h = data->h;

    double **u_new = make2DDoubleArray(m, n);
    double **v_new = make2DDoubleArray(m, n);

    double **u_n_1 = data->u_n_1;
    double **u_n = data->u_n;
    double **v_n_1 = data->v_n_1;
    double **v_n = data->v_n;
    double **P = data->P;

    double **u_star = data->u_star;
    double **v_star = data->v_star;
    double **phi = data->phi;

    int i, j;
    /* Set boundary for u_new, v_new */

    for (j = 1; j < n - 1; j++) {
        u_new[0][j] = u_n[0][j];
    }
    for (i = 1; i < m - 1; i++) {
        for (j = 1; j < n - 1; j++) {
            u_new[i][j] = u_star[i][j] - dt * (phi[i + 1][j] - phi[i][j]) / h;
            v_new[i][j] = v_star[i][j] - dt * (phi[i][j + 1] - phi[i][j]) / h;
            P[i][j] += phi[i][j];
        }
    }

    free2Darray(u_n_1, m);
    free2Darray(v_n_1, m);

    data->u_n_1 = u_n;
    data->u_n = u_new;

    data->v_n_1 = v_n;
    data->v_n = v_new;

}

void update_scalars(Data* data)
{
    /* TEMPERATURE AND SPECIES */

    int m = data->m;
    int n = data->n;
    double h = data->h;
    int i,j;

    double **T_new = make2DDoubleArray(m, n);
    double **T_n = data->T_n;
    double **T_n_1 = data->T_n_1;
    double **Ts = data->Ts;

    double **u_n = data->u_n;
    double **v_n = data->v_n;

    int Ns = data->Ns;
    double ***C_new = make3DDoubleArray(Ns, m, n);
    double ***C_n = data->C_n;
    double ***C_n_1 = data->C_n_1;
    double ***Cs = data->Cs;

    double alpha_f = data->alpha_f;
    double *Df = data->Df;
    double **I_S = data->I_S;
    double ramp = data->ramp;
    double dtau = data->dtau;
    double dt = data->dt;

    double H_T, H_T_old, lapT;
    double H_C, H_C_old, lapC;

    /** TEMPERATURE **/

    for (i = 1; i < m - 1; i++) {
        for (j = 1; j < n - 1; j++) {

            // ADVECTIVE TERMS
            H_T = .5*(u_n[i][j]*(T_n[i+1][j]-T_n[i][j])/h + u_n[i-1][j]*(T_n[i][j]-T_n[i-1][j])/h)
                  + .5*(v_n[i][j]*(T_n[i][j+1]-T_n[i][j])/h + v_n[i][j-1]*(T_n[i][j]-T_n[i][j-1])/h);

            if(data->iter == 1)
            {
                H_T_old = H_T;
            }
            else{
                H_T_old = data->H_T_n_1[i][j];
            }

            // DIFFUSION TERM
            lapT = (T_n[i + 1][j] + T_n[i - 1][j] + T_n[i][j + 1] + T_n[i][j - 1] - 4. * T_n[i][j]) / (h * h);

#ifdef EXPLICIT
            // EXPLICIT VERSION
            T_new[i][j] = T_n[i][j] + dt * (-1.5 * H_T + 0.5 * H_T_old
                                            + alpha_f * lapT)
                                    - ramp*I_S[i][j]*(dt/dtau)*(T_n[i][j] - Ts[i][j]);
#else
            // IMPLICIT VERSION
            T_new[i][j] = (T_n[i][j] + dt * (-1.5 * H_T + 0.5 * H_T_old
                                             + alpha_f * lapT
                                             + ramp * I_S[i][j] * Ts[i][j] / dtau)) /
                          (1. + ramp * I_S[i][j] * dt / dtau);
#endif
            data->H_T_n_1[i][j] = H_T;

            /** SPECIES **/

            for (int s = 0; s < Ns-1; s++) {
                // ADVECTIVE TERMS
                H_C = .5*(u_n[i][j]*(C_n[s][i+1][j]-C_n[s][i][j])/h + u_n[i-1][j]*(C_n[s][i][j]-C_n[s][i-1][j])/h)
                      + .5*(v_n[i][j]*(C_n[s][i][j+1]-C_n[s][i][j])/h + v_n[i][j-1]*(C_n[s][i][j]-C_n[s][i][j-1])/h);

                if(data->iter == 1)
                {
                    H_C_old = H_C;
                }
                else{
                    H_C_old = data->H_C_n_1[s][i][j];
                }

                // DIFFUSION TERM
                lapC = (C_n[s][i + 1][j] + C_n[s][i - 1][j] + C_n[s][i][j + 1] + C_n[s][i][j - 1] - 4. * C_n[s][i][j]) / (h * h);
#ifdef EXPLICIT
                //EXPLICIT VERSION
              C_new[s][i][j] = C_n[s][i][j] + dt * (-1.5 * H_C + 0.5 * H_C_old
                                                     + Df[s] * lapC)
                                              - ramp*I_S[i][j]*(dt/dtau)*(C_n[s][i][j] - Cs[s][i][j]);
#else
                // IMPLICIT VERSION
                C_new[s][i][j] = (C_n[s][i][j] + dt * (-1.5 * H_C + 0.5 * H_C_old
                                                       + Df[s] * lapC
                                                       + ramp * I_S[i][j] * Cs[s][i][j] / dtau))/
                                 (1. + ramp * I_S[i][j] * dt / dtau);

#endif
                data->H_C_n_1[s][i][j] = H_C;
            }

        }
    }

    free2Darray(T_n_1, m);
    free3Darray(C_n_1, Ns, m);

    data->T_n_1 = T_n;
    data->T_n = T_new;

    data->C_n_1 = C_n;
    data->C_n = C_new;
}

void get_vorticity(Data* data){
    int m = data->m;
    int n = data->n;
    double h = data->h;

    int i,j;

    double** omega = data->omega;
    double** u_n = data->u_n;
    double** v_n = data->v_n;

    for (i = 0; i < m - 1; i++) {
        for (j = 1; j < n - 2; j++) {
            omega[i][j] = (v_n[i + 1][j] - v_n[i][j]) / h - (u_n[i][j + 1] - u_n[i][j]) / h;
        }
    }
}

