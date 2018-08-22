//
// Created by Baptiste Hardy on 5/07/18.
//

#include "main.h"
#include "flow_solver.h"

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
                xG = xg[k][2];
                yG = yg[k][2];

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

                if((int) floor((delta-theta[k][2])/(M_PI/2.)) % 2 == 0 ){
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
            Cs[0][i][j] = 0.;
            for(int k=0; k<Np; k++){
                Cs[0][i][j] += chi_S[k][i][j]*Cp[k][0];
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
                u_s[i][j]+= chi_U[k][i][j]*(Up[k][2] - Omega_p[k][2]*(yU-yg[k][2]));
                v_s[i][j]+= chi_V[k][i][j]*(Vp[k][2] + Omega_p[k][2]*(xV-xg[k][2]));
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

    double** T_n = data->T_n;

    double Um = data->u_m;
    double Gr = data->Gr;
    double Re_p = data->Re_p;

    double H_U, H_U_old;
    double H_V, H_V_old;
    double lapU, lapV;
    double dpdx, dpdy;

    double T_u = 0;
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

#ifdef TEMP
            //BUOYANCY TERM
            T_u = .5*(T_n[i][j] + T_n[i+1][j]);
#endif

#ifdef EXPLICIT
            // EXPLICIT VERSION
            u_star[i][j] = u_n[i][j] + dt*(-1.5*H_U + 0.5*H_U_old - dpdx + nu*lapU + Gr/pow(Re_p,2)*T_n[i][j]) - ramp*I_U[i][j]*(dt/dtau)*(u_n[i][j] - u_s[i][j]);
#else
            // IMPLICIT VERSION
            u_star[i][j] = (u_n[i][j] + dt*(-1.5*H_U + 0.5*H_U_old - dpdx + nu*lapU + (Gr/pow(Re_p,2))*T_u) + (dt/dtau)*ramp*I_U[i][j]*u_s[i][j])/(1.+ramp*I_U[i][j]*dt/dtau);
#endif
            data->H_u_n_1[i][j] = H_U;
        }
    }

    /*Outflow condition */
//    for (j=1; j<n-1; j++){
//        u_star[m-2][j] = u_n[m-2][j] - dt*Um*(u_n[m-2][j]-u_n[m-3][j])/h;
//    }
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
    int Np = data->Np;
    double h = data->h;
    double dx = data->h;
    double dy = data->h;
    int i,j, k;

    double ramp = data->ramp;
    double dtau = data->dtau;
    double dt = data->dt;

    double rho_f = data->rho_f;
    double rho_s = data->rho_s;
    double cp_s = data->cp_s;
    double cp_f = data->cp_f;
    double dH = data->dH;

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
    double** kappa = data->kappa;
    double *Df = data->Df;
    double ***D = data->D;
    double rho_star, cp_star;

    double **I_S = data->I_S;
    double q_left, q_right, q_top, q_bottom;
    double j_left, j_right, j_top, j_bottom;

    double H_T_n, H_T_n_1, diff_T, S_T;
    double H_C_n, H_C_n_1, diff_C, S_C;

    double* rate = make1DDoubleArray(Ns);

    /** TEMPERATURE **/

    for (i = 1; i < m - 1; i++) {
        for (j = 1; j < n - 1; j++) {


            get_rate(data, rate, C_n, T_n, i, j);

            // ADVECTIVE TERMS
            H_T_n = .5*(u_n[i][j]*(T_n[i+1][j]-T_n[i][j])/h + u_n[i-1][j]*(T_n[i][j]-T_n[i-1][j])/h)
                  + .5*(v_n[i][j]*(T_n[i][j+1]-T_n[i][j])/h + v_n[i][j-1]*(T_n[i][j]-T_n[i][j-1])/h);

            if(data->iter == 1)
            {
                H_T_n_1 = H_T_n;
            }
            else{
                H_T_n_1 = data->H_T_n_1[i][j];
            }

            // DIFFUSION TERM
#ifdef INTRAPARTICLE

            int K[4] = {-1, -1, -1, -1};
            double THETA[4] = {0., 0., 0., 0.};
            double kappa_right, kappa_left, kappa_top, kappa_bottom;
            track_interface(data, K, THETA, i, j);

            kappa_right = kappa[i][j]*kappa[i+1][j]/(kappa[i+1][j]*THETA[0] + kappa[i][j]*(1.-THETA[0]));
            kappa_left = kappa[i][j]*kappa[i-1][j]/(kappa[i-1][j]*THETA[1] + kappa[i][j]*(1.-THETA[1]));
            kappa_top = kappa[i][j]*kappa[i][j+1]/(kappa[i][j+1]*THETA[2] + kappa[i][j]*(1.-THETA[2]));
            kappa_bottom = kappa[i][j]*kappa[i][j-1]/(kappa[i][j-1]*THETA[3] + kappa[i][j]*(1.-THETA[3]));

            q_right = -kappa_right*(T_n[i+1][j] - T_n[i][j])/dx;
            q_left = -kappa_left*(T_n[i][j] - T_n[i-1][j])/dx;
            q_top = -kappa_top*(T_n[i][j+1] - T_n[i][j])/dy;
            q_bottom = -kappa_bottom*(T_n[i][j] - T_n[i][j-1])/dy;

            rho_star = 1; // (1-I_S[i][j])*rho_f + I_S[i][j]*rho_s;
            cp_star = 1; //(1-I_S[i][j])*cp_f + I_S[i][j]*cp_s;

            diff_T = -((q_right-q_left)/dx + (q_top-q_bottom)/dy)/(rho_star*cp_star);

            S_T = fabs(rate[0])*(-dH)/(rho_star*cp_star);
            T_new[i][j] = T_n[i][j] + dt*((-1.5*H_T_n_1 + 0.5*H_T_n_1) + diff_T + I_S[i][j]*S_T);
#else
            diff_T = (T_n[i + 1][j] + T_n[i - 1][j] + T_n[i][j + 1] + T_n[i][j - 1] - 4. * T_n[i][j]) / (h * h);

#ifdef EXPLICIT
            // EXPLICIT VERSION
            T_new[i][j] = T_n[i][j] + dt * (-1.5 * H_T_n + 0.5 * H_T_n_1
                                            + alpha_f * diff_T)
                                    - ramp*I_S[i][j]*(dt/dtau)*(T_n[i][j] - Ts[i][j]);
#else
            // IMPLICIT VERSION
            T_new[i][j] = (T_n[i][j] + dt * (-1.5 * H_T_n + 0.5 * H_T_n_1
                                             + alpha_f * diff_T
                                             + ramp * I_S[i][j] * Ts[i][j] / dtau)) /
                          (1. + ramp * I_S[i][j] * dt / dtau);
#endif

#endif
            data->H_T_n_1[i][j] = H_T_n;

            /** SPECIES **/

            for (int s = 0; s < Ns; s++) {

                // ADVECTIVE TERMS
                H_C_n = .5*(u_n[i][j]*(C_n[s][i+1][j]-C_n[s][i][j])/h + u_n[i-1][j]*(C_n[s][i][j]-C_n[s][i-1][j])/h)
                      + .5*(v_n[i][j]*(C_n[s][i][j+1]-C_n[s][i][j])/h + v_n[i][j-1]*(C_n[s][i][j]-C_n[s][i][j-1])/h);

                if(data->iter == 1)
                {
                    H_C_n_1 = H_C_n;
                }
                else{
                    H_C_n_1 = data->H_C_n_1[s][i][j];
                }

                // DIFFUSION TERM

#ifdef INTRAPARTICLE

                j_left = -.5*(D[s][i][j] + D[s][i-1][j])*(C_n[s][i][j] - C_n[s][i-1][j])/dx;
                j_right = -.5*(D[s][i+1][j] + D[s][i][j])*(C_n[s][i+1][j] - C_n[s][i][j])/dx;
                j_top = -.5*(D[s][i][j+1] + D[s][i][j])*(C_n[s][i][j+1] - C_n[s][i][j])/dy;
                j_bottom = -.5*(D[s][i][j] + D[s][i][j-1])*(C_n[s][i][j] - C_n[s][i][j-1])/dy;

                diff_C = -( (j_right-j_left)/dx + (j_top-j_bottom)/dy );

                C_new[s][i][j] = C_n[s][i][j] + dt*((-1.5*H_C_n + 0.5*H_C_n_1) + diff_C + I_S[i][j]*rate[s]);

#else
                diff_C = (C_n[s][i + 1][j] + C_n[s][i - 1][j] + C_n[s][i][j + 1] + C_n[s][i][j - 1] - 4. * C_n[s][i][j]) / (h * h);
#ifdef EXPLICIT
                //EXPLICIT VERSION
              C_new[s][i][j] = C_n[s][i][j] + dt * (-1.5 * H_C_n + 0.5 * H_C_n_1
                                                     + Df[s] * diff_C)
                                              - ramp*I_S[i][j]*(dt/dtau)*(C_n[s][i][j] - Cs[s][i][j]);
#else
                // IMPLICIT VERSION
                C_new[s][i][j] = (C_n[s][i][j] + dt * (-1.5 * H_C_n + 0.5 * H_C_n_1
                                                       + Df[s] * diff_C
                                                       + ramp * I_S[i][j] * Cs[s][i][j] / dtau))/
                                 (1. + ramp * I_S[i][j] * dt / dtau);

#endif
#endif
                data->H_C_n_1[s][i][j] = H_C_n;
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

void get_vorticity(Data* data)
{
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

void get_rate(Data* data, double* r, double*** C, double** Ts, int i, int j)
{
    double k = data->Da;
    r[0] = -k*C[0][i][j];
    //r[1] = -r[0];
}

void get_diffusivity(Data* data)
{
    double*** D = data->D;
    double** I_S = data->I_S;
    double* Df = data->Df;
    double Ds;

    int i, j, s;
    int m = data->m;
    int n = data->n;
    int Ns = data->Ns;

    for(i=0; i<m; i++) {
        for(j=0; j<n; j++){
            for(s = 0; s<Ns; s++) {
                Ds = 0.1*Df[s];
                D[s][i][j] = (1 - I_S[i][j])*Df[s] + I_S[i][j] * Ds;
            }
        }
    }
}

void get_conductivity(Data* data)
{
    double** I_S = data->I_S;
    double** kappa = data->kappa;
    double kappa_f = data->kappa_f;
    double kappa_s = kappa_f;

    int i, j;
    int m = data->m;
    int n = data->n;

    for(i=0; i<m; i++) {
        for(j=0; j<n; j++){
            kappa[i][j] = (1-I_S[i][j])*kappa_f + I_S[i][j]*kappa_s;
        }
    }
}

void track_interface(Data* data, int* K, double* THETA, int i, int j)
{
    int k;
    int Np = data->Np;
    double h = data->h;

    double** xg = data->xg;
    double** yg = data->yg;
    double* rp = data->rp;
    double** I_S = data->I_S;
    double*** Ip_S = data->Ip_S;

    double xi, yj, X_I_1, X_I_2, Y_I_1, Y_I_2, dx1, dx2, dy1, dy2, d;

    yj = (j-0.5) * h;
    xi = (i-0.5) * h;

    /** RIGHT ARM **/
    if (I_S[i][j] == I_S[i+1][j])
    {
        // NO INTERFACE BETWEEN i and i+1
        THETA[0] = 0.5;
    }
    else {
        for (k = 0; k < Np; k++) {
            if (Ip_S[k][i][j] != Ip_S[k][i+1][j])
            {
                X_I_1 = xg[0][2] - sqrt(pow(rp[k], 2) - pow(yj - yg[k][2], 2));
                X_I_2 = xg[0][2] + sqrt(pow(rp[k], 2) - pow(yj - yg[k][2], 2));

                dx1 = X_I_1 - xi;
                dx2 = X_I_2 - xi;

                d = fmin(fabs(dx1), fabs(dx2));
                THETA[0] = d/h;

                break;
            }
        }
    }

    /** LEFT ARM **/
    if (I_S[i][j] == I_S[i-1][j])
    {
        // NO INTERFACE BETWEEN i and i-1
        THETA[1] = 0.5;
    }
    else {
        for (k = 0; k < Np; k++) {
            if (Ip_S[k][i][j] != Ip_S[k][i-1][j])
            {
                X_I_1 = xg[0][2] - sqrt(pow(rp[k], 2) - pow(yj - yg[k][2], 2));
                X_I_2 = xg[0][2] + sqrt(pow(rp[k], 2) - pow(yj - yg[k][2], 2));

                dx1 = X_I_1 - xi;
                dx2 = X_I_2 - xi;

                d = fmin(fabs(dx1), fabs(dx2));
                THETA[1] = d/h;

                break;
            }
        }
    }

    /** TOP ARM **/
    if (I_S[i][j] == I_S[i][j+1])
    {
        // NO INTERFACE BETWEEN j and j+1
        THETA[2] = 0.5;
    }
    else {
        for (k = 0; k < Np; k++) {
            if (Ip_S[k][i][j] != Ip_S[k][i][j+1])
            {
                Y_I_1 = yg[0][2] - sqrt(pow(rp[k], 2) - pow(xi - xg[k][2], 2));
                Y_I_2 = yg[0][2] + sqrt(pow(rp[k], 2) - pow(xi - xg[k][2], 2));

                dy1 = Y_I_1 - yj;
                dy2 = Y_I_2 - yj;

                d = fmin(fabs(dy1), fabs(dy2));
                THETA[2] = d/h;
                break;
            }
        }
    }

    /** BOTTOM ARM **/
    if (I_S[i][j] == I_S[i][j-1])
    {
        // NO INTERFACE BETWEEN j and j-1
        THETA[3] = 0.5;
    }
    else {
        for (k = 0; k < Np; k++) {
            if (Ip_S[k][i][j] != Ip_S[k][i][j-1])
            {
                Y_I_1 = yg[0][2] - sqrt(pow(rp[k], 2) - pow(xi - xg[k][2], 2));
                Y_I_2 = yg[0][2] + sqrt(pow(rp[k], 2) - pow(xi - xg[k][2], 2));

                dy1 = Y_I_1 - yj;
                dy2 = Y_I_2 - yj;

                d = fmin(fabs(dy1), fabs(dy2));
                THETA[3] = d/h;

                break;
            }
        }
    }


//
//        if ( (dx1 >= 0 && dx1 < h) || (dx2 >= 0 && dx2 < h) ){
//            *right = 1;
//            K[0] = k;
//            THETA[0] = fmin(dx1, dx2)/h;
//        }
//        if ( (dx1 < 0 && dx1 >= -h) || (dx2 < 0 && dx2 >= -h) ){
//            *left = 1;
//            K[1] = k;
//            THETA[1] = fmin(fabs(dx1), fabs(dx2))/h;
//        }
//        if ( (dy1 >= 0 && dy1 < h) || (dy2 >= 0 && dy2 < h) ){
//            *above = 1;
//            K[2] = k;
//            THETA[2] = fmin(dy1, dy2)/h;
//        }
//        if ( (dy1 < 0 && dy1 >= -h) || (dy2 < 0 && dy2 >= -h) ){
//            *below = 1;
//            K[3] = k;
//            THETA[3] = fmin(fabs(dy1), fabs(dy2))/h;
//        }

}



double get_tg_gradient(Data* data, int i, int j, int k)
{
    double X = data->xg[k][2];
    double Y = data->yg[k][2];


}