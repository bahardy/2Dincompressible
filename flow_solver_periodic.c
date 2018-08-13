//
// Created by Baptiste Hardy on 5/07/18.
//

#include "flow_solver_periodic.h"


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
        /* Bottom and top Wall : slip : du/dn= 0 */
        u_n[i][0] = u_n[i][1];
        u_n[i][n-1] = u_n[i][n-2];

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
    double* rp = data->rp;
    double dist;
    int m = data->m;
    int n = data->n;
    double h = data->h;
    double L = data->L;
    int Np = data->Np;
    int b1, b2, b3;
    double d1, d2, d3;
    double xU, xV, yU, yV, yS, xS;
    double** xg = data->xg;
    double** yg = data->yg;
    double** theta = data->theta;

    double xloc, yloc, delta;


    for(int i=0; i<m; i++){
        xU = (i+1)*h;
        xV = (i+0.5)*h;
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
                double xG = fmod(xg[k][2], L);
                double yG = yg[k][2];

#ifdef ELLIPSE
                double x;
		double y;
		/*ELLIPSE*/
		x = xU-Xp_k[k];
		y = yU-yg[k];
                double EU = pow(b,2)*(pow(y*cos(theta[k]),2)+pow(x*sin(theta[k]),2)-y*x*sin(2*theta[k]))
                + pow(a,2)*(pow(x*cos(theta[k]),2)+pow(y*sin(theta[k]),2)+y*x*sin(2*theta[k]))
                - pow(a,2)*pow(b,2);

		x = xV-Xp_k[k];
		y = yV-yG;
                double EV = pow(b,2)*(pow(y*cos(theta[k]),2)+pow(x*sin(theta[k]),2)-y*x*sin(2*theta[k]))
                + pow(a,2)*(pow(x*cos(theta[k]),2)+pow(y*sin(theta[k]),2)+y*x*sin(2*theta[k]))
                - pow(a,2)*pow(b,2);

		x = xS-Xp_k[k];
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
                b1 = ((xS-xG)*(xS-xG)+(yS-yG)*(yS-yG)<= (rp[k]+data->eps)*(rp[k]+data->eps));
                b2 = ((xS-(xG+L))*(xS-(xG+L))+(yS-yG)*(yS-yG)<= (rp[k]+data->eps)*(rp[k]+data->eps));
                b3 = ((xS-(xG-L))*(xS-(xG-L))+(yS-yG)*(yS-yG)<= (rp[k]+data->eps)*(rp[k]+data->eps));
                chi_S[k][i][j]= (b1 || b2 || b3);

                b1 = ((xU-xG)*(xU-xG)+(yU-yG)*(yU-yG)<= (rp[k]+data->eps)*(rp[k]+data->eps));
                b2 = ((xU-(xG+L))*(xU-(xG+L))+(yU-yG)*(yU-yG)<= (rp[k]+data->eps)*(rp[k]+data->eps));
                b3 = ((xU-(xG-L))*(xU-(xG-L))+(yU-yG)*(yU-yG)<= (rp[k]+data->eps)*(rp[k]+data->eps));

                chi_U[k][i][j]= (b1 || b2 || b3);

                b1 = ((xV-xG)*(xV-xG)+(yV-yG)*(yV-yG)<= (rp[k]+data->eps)*(rp[k]+data->eps));
                b2 = ((xV-(xG+L))*(xV-(xG+L))+(yV-yG)*(yV-yG)<= (rp[k]+data->eps)*(rp[k]+data->eps));
                b3 = ((xV-(xG-L))*(xV-(xG-L))+(yV-yG)*(yV-yG)<= (rp[k]+data->eps)*(rp[k]+data->eps));

                chi_V[k][i][j]= (b1 || b2 || b3);

#ifndef SMOOTHING
                b1 = ((xS-xG)*(xS-xG)+(yS-yG)*(yS-yG)<= rp[k]*rp[k]);
                b2 = ((xS-(xG+L))*(xS-(xG+L))+(yS-yG)*(yS-yG)<= rp[k]*rp[k]);
                Ip_S[k][i][j]= b1 || b2;

                b1 = ((xU-xG)*(xU-xG)+(yU-yG)*(yU-yG)<= rp[k]*rp[k]);
                b2 = ((xU-(xG+L))*(xU-(xG+L))+(yU-yG)*(yU-yG)<= rp[k]*rp[k]);
                Ip_U[k][i][j]= b1 || b2;

                b1 = ((xV-xG)*(xV-xG)+(yV-yG)*(yV-yG)<= rp[k]*rp[k]);
                b2 = ((xV-(xG+L))*(xV-(xG+L))+(yV-yG)*(yV-yG)<= rp[k]*rp[k]);
                Ip_V[k][i][j]= b1 || b2;

#endif

#ifdef SMOOTHING

                //Smoothing S
                d1 = sqrt((xS-xG)*(xS-xG)+(yS-yG)*(yS-yG));
                d2 = sqrt((xS-(xG+L))*(xS-(xG+L))+(yS-yG)*(yS-yG));
                d3 = sqrt((xS-(xG-L))*(xS-(xG-L))+(yS-yG)*(yS-yG));
                dist =  rp[k] - fmin(d1, fmin(d2,d3));

                if( dist < - data->eps)
                    Ip_S[k][i][j] = 0;
                else if( fabs(dist) <= data->eps)
                    Ip_S[k][i][j] = .5*(1 + dist/data->eps + (1./M_PI)*sin( M_PI* dist/data->eps) );
                else if( dist > data->eps)
                    Ip_S[k][i][j] = 1;

                //Smoothing U
                d1 = sqrt((xU-xG)*(xU-xG)+(yU-yG)*(yU-yG));
                d2 = sqrt((xU-(xG+L))*(xU-(xG+L))+(yU-yG)*(yU-yG));
                d3 = sqrt((xU-(xG-L))*(xU-(xG-L))+(yU-yG)*(yU-yG));
                dist =  rp[k] - fmin(d1, fmin(d2,d3));

                if( dist < - data->eps)
                    Ip_U[k][i][j] = 0;
                else if( fabs(dist) <=data->eps)
                    Ip_U[k][i][j] = .5*(1 + dist/data->eps + (1./M_PI)*sin( M_PI* dist/data->eps) );
                else if( dist > data->eps)
                    Ip_U[k][i][j] = 1;


                //Smoothing V
                d1 = sqrt((xV-xG)*(xV-xG)+(yV-yG)*(yV-yG));
                d2 = sqrt((xV-(xG+L))*(xV-(xG+L))+(yV-yG)*(yV-yG));
                d3 = sqrt((xV-(xG-L))*(xV-(xG-L))+(yV-yG)*(yV-yG));
                dist =  rp[k] - fmin(d1, fmin(d2,d3));

                if( dist < - data->eps)
                    Ip_V[k][i][j] = 0;
                else if( fabs(dist) <= data->eps)
                    Ip_V[k][i][j] = .5*(1 + dist/data->eps + (1./M_PI)*sin( M_PI* dist/data->eps) );
                else if( dist > data->eps)
                    Ip_V[k][i][j] = 1;

#endif

                xloc = xS-xg[k][2];
                yloc = yS-yG;
                delta = atan2(yloc, xloc);
                coloring[i][j] += Ip_S[k][i][j];

                if((int) floor(((delta-theta[k][1])/(M_PI/2.))) % 2 == 0 ){
                    coloring[i][j] = -coloring[i][j];
                }
#endif
                I_S[i][j] += Ip_S[k][i][j];
                I_U[i][j] += Ip_U[k][i][j];
                I_V[i][j] += Ip_V[k][i][j];
            }

            if(I_S[i][j] > 1 || I_U[i][j] > 1 || I_V[i][j] > 1 ){
                PetscPrintf(PETSC_COMM_WORLD, "Collision de particules \n");
                I_S[i][j] = fmin(I_S[i][j], 1);
                I_U[i][j] = fmin(I_U[i][j], 1);
                I_V[i][j] = fmin(I_V[i][j], 1);
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
    double L = data->L;
    int Np = data->Np;
    double h = data->h;

    double xV, yU;
    for(int i=0; i<m; i++){
        xV = (i+0.5)*h;
        for(int j=0; j<n; j++){
            yU = (j-0.5)*h;
            u_s[i][j] = 0.;
            v_s[i][j] = 0.;
            for (int k = 0; k<Np; k++){
                u_s[i][j]+= chi_U[k][i][j]*(Up[k][2] - Omega_p[k][2]*(yU-yg[k][2]) );
                v_s[i][j]+= chi_V[k][i][j]*(Vp[k][2] + Omega_p[k][2]*fmod(xV-xg[k][2],L) );
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
    double** u_n_1 = data->u_n_1;
    double** u_star = data->u_star;
    double** u_s = data->u_s;

    double** v_n = data->v_n;
    double** v_n_1 = data->v_n_1;
    double** v_star = data->v_star;
    double** v_s = data->v_s;

    double** P = data->P;

    double H_U, H_U_old;
    double H_V, H_V_old;
    double lapU, lapV;
    double dpdx, dpdy;

    double uL, uR, vT, vB;
    double dudxR, dudxL, dudyT, dudyB, dvdxR, dvdxL, dvdyT, dvdyB;

    /* u_star ADAMS-BASHFORTH 2 */
    for (i=0; i<m; i++){
        for (j=1; j<n-1; j++){

            // CONVECTIVE TERM
            /** time n-1 **/
            uR = .5*(u_n_1[i][j] + u_n_1[(i+1+m)%m][j]);
            uL = .5*(u_n_1[(i-1+m)%m][j] + u_n_1[i][j]);
            dudxR = (u_n_1[(i+1+m)%m][j]- u_n_1[i][j])/h;
            dudxL = (u_n_1[i][j]- u_n_1[(i-1+m)%m][j])/h;

            vT = .5*(v_n_1[i][j] + v_n_1[(i+1+m)%m][j]);
            vB = .5*(v_n_1[i][j-1] + v_n_1[(i+1+m)%m][j-1]);
            dudyT = (u_n_1[i][j+1] - u_n_1[i][j])/h;
            dudyB = (u_n_1[i][j] - u_n_1[i][j-1])/h;

            H_U_old = .5*(uR*dudxR + uL*dudxL) + .5*(vT*dudyT + vB*dudyB);

            /** time n **/

            uR = .5*(u_n[i][j] + u_n[(i+1+m)%m][j]);
            uL = .5*(u_n[(i-1+m)%m][j] + u_n[i][j]);
            dudxR = (u_n[(i+1+m)%m][j]- u_n[i][j])/h;
            dudxL = (u_n[i][j]- u_n[(i-1+m)%m][j])/h;

            vT = .5*(v_n[i][j] + v_n[(i+1+m)%m][j]);
            vB = .5*(v_n[i][j-1] + v_n[(i+1+m)%m][j-1]);
            dudyT = (u_n[i][j+1] - u_n[i][j])/h;
            dudyB = (u_n[i][j] - u_n[i][j-1])/h;

            H_U = .5*(uR*dudxR + uL*dudxL) + .5*(vT*dudyT + vB*dudyB);


            // LAPLACIAN
            lapU = (u_n[(i+1+m)%m][j]+u_n[(i-1+m)%m][j]+u_n[i][j+1]+u_n[i][j-1]-4.*u_n[i][j])/(h*h);

            // PRESSURE TERM
            dpdx = (P[(i+1+m)%m][j]-P[i][j])/h;

#ifdef EXPLICIT
            //EXPLICIT VERSION
            u_star[i][j] = u_n[i][j] + dt*(-1.5*H_U + 0.5*H_U_old - dpdx + nu*lapU - ramp*I_U[i][j]*(u_n[i][j] - u_s[i][j])/dtau);

#else
            //IMPLICIT VERSION
            u_star[i][j] = (u_n[i][j] + dt*(-1.5*H_U + 0.5*H_U_old - dpdx + nu*lapU) + (dt/dtau)*ramp*I_U[i][j]*u_s[i][j])/(1.+ramp*I_U[i][j]*dt/dtau);
#endif

        }
    }

    /* v_star  ADAMS-BASHFORTH 2 */
    for (i=0; i<m; i++){
        for (j=1; j<n-2; j++){

            // CONVECTIVE TERM

            uR = .5*(u_n_1[i][j] + u_n_1[i][j+1]);
            uL = .5*(u_n_1[(i-1+m)%m][j] + u_n_1[(i-1+m)%m][j+1]);
            dvdxR = (v_n_1[(i+1+m)%m][j]- v_n_1[i][j])/h;
            dvdxL = (v_n_1[i][j]- v_n_1[(i-1+m)%m][j])/h;

            vT = .5*(v_n_1[i][j] + v_n_1[i][j+1]);
            vB = .5*(v_n_1[i][j] + v_n_1[i][j-1]);
            dvdyT = (v_n_1[i][j+1] - v_n_1[i][j])/h;
            dvdyB = (v_n_1[i][j] - v_n_1[i][j-1])/h;

            H_V_old = .5*(uR*dvdxR + uL*dvdxL) + .5*(vT*dvdyT + vB*dvdyB);

            uR = .5*(u_n[i][j] + u_n[i][j+1]);
            uL = .5*(u_n[(i-1+m)%m][j] + u_n[(i-1+m)%m][j+1]);
            dvdxR = (v_n[(i+1+m)%m][j]- v_n[i][j])/h;
            dvdxL = (v_n[i][j]- v_n[(i-1+m)%m][j])/h;

            vT = .5*(v_n[i][j] + v_n[i][j+1]);
            vB = .5*(v_n[i][j] + v_n[i][j-1]);
            dvdyT = (v_n[i][j+1] - v_n[i][j])/h;
            dvdyB = (v_n[i][j] - v_n[i][j-1])/h;

            H_V = .5*(uR*dvdxR + uL*dvdxL) + .5*(vT*dvdyT + vB*dvdyB);

            // LAPLACIAN
            lapV = (v_n[(i+1+m)%m][j]+v_n[(i-1+m)%m][j]+v_n[i][j+1]+v_n[i][j-1]-4.*v_n[i][j])/(h*h);

            // PRESSURE TERM
            dpdy = (P[i][j+1]-P[i][j])/h;


#ifdef EXPLICIT
            //EXPLICIT VERSION
            v_star[i][j] = v_n[i][j] + dt*(-1.5*H_V + 0.5*H_V_old - dpdy + nu*lapV - ramp*I_V[i][j]*(v_n[i][j] - v_s[i][j])/dtau);
#else
            //IMPLICIT VERSION
            v_star[i][j] = (v_n[i][j] + dt*(-1.5*H_V + 0.5*H_V_old - dpdy + nu*lapV) + (dt/dtau)*ramp*I_V[i][j]*v_s[i][j])/(1.+ramp*I_V[i][j]*dt/dtau);
#endif
        }
    }

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

    for (i = 0; i < m; i++) {
        for (j = 1; j < n - 2; j++) {
            omega[i][j] = (v_n[(i+1+m)%m][j] - v_n[i][j]) / h - (u_n[i][j + 1] - u_n[i][j]) / h;
        }
    }
}

void update_flow(Data* data)
{

    int m = data->m;
    int n = data->n;
    double dt = data->dt;
    double h = data->h;
    //double H = data->H;

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

    for (i = 0; i < m; i++) {
        for (j = 0; j < n - 1; j++) {
            u_new[i][j] = u_star[i][j] - dt * (phi[(i+1 + m) % m][j] - phi[i][j]) / h;
            v_new[i][j] = v_star[i][j] - dt * (phi[i][j+1] - phi[i][j]) / h;
            P[i][j] += phi[i][j];
        }
    }

#ifdef TEMP


#endif

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
    int i,j;
    int m = data->m;
    int n = data->n;
    double dt = data->dt;
    double h = data->h;

    double **u_n_1 = data->u_n_1;
    double **u_n = data->u_n;
    double **v_n_1 = data->v_n_1;
    double **v_n = data->v_n;

    double **T_new = make2DDoubleArray(m, n);
    double **T_n = data->T_n;
    double **T_n_1 = data->T_n_1;
    double **Ts = data->Ts;
    double ***Cs = data->Cs;

    int Ns = data->Ns;
    double ***C_new = make3DDoubleArray(Ns, m, n);
    double ***C_n = data->C_n;
    double ***C_n_1 = data->C_n_1;

    double alpha_f = data->alpha_f;
    double *Df = data->Df;
    double **I_S = data->I_S;
    double ramp = data->ramp;
    double dtau = data->dtau;

    double Uij, Uij_old, Vij, Vij_old;
    double H_T, H_T_old, lapT;
    double H_C, H_C_old, lapC;

    for (i = 1; i < m - 1; i++) {
        for (j = 1; j < n - 1; j++) {

            // Need to have a value for velocities at cell center
            Uij = .5 * (u_n[i][j] + u_n[i - 1][j]); /* average of the two neighbour points along x-direction */
            Vij = .5 * (v_n[i][j] + v_n[i - 1][j]); /* average of the two neighbour points along y-direction */
            Uij_old = .5 * (u_n_1[i][j] + u_n_1[i - 1][j]);
            Vij_old = .5 * (v_n_1[i][j] + v_n_1[i - 1][j]);

            // Advective terms
            H_T_old = Uij_old * (T_n_1[i + 1][j] - T_n_1[i - 1][j]) / (2. * h) +
                      Vij_old * (T_n_1[i][j + 1] - T_n_1[i][j - 1]) / (2. * h);
            H_T = Uij * (T_n[i + 1][j] - T_n[i - 1][j]) / (2. * h) + Vij * (T_n[i][j + 1] - T_n[i][j - 1]) / (2. * h);
            // Laplacian
            lapT = (T_n[i + 1][j] + T_n[i - 1][j] + T_n[i][j + 1] + T_n[i][j - 1] - 4. * T_n[i][j]) / (h * h);

            T_new[i][j] = (T_n[i][j] + dt * (-1.5 * H_T + 0.5 * H_T_old
                                             + alpha_f * lapT
                                             + ramp * I_S[i][j] * Ts[i][j] / dtau)) /
                          (1. + ramp * I_S[i][j] * dt / dtau);

            for (int s = 0; s < Ns; s++) {
                // Advective terms
                H_C_old = Uij_old * (C_n_1[s][i + 1][j] - C_n_1[s][i - 1][j]) / (2. * h) +
                          Vij_old * (C_n_1[s][i][j + 1] - C_n_1[s][i][j - 1]) / (2. * h);
                H_C = Uij * (C_n[s][i + 1][j] - C_n[s][i - 1][j]) / (2. * h) +
                      Vij * (C_n[s][i][j + 1] - C_n[s][i][j - 1]) / (2. * h);
                // Laplacian
                lapC = (C_n[s][i + 1][j] + C_n[s][i - 1][j] + C_n[s][i][j + 1] + C_n[s][i][j - 1] - 4. * C_n[s][i][j]) / (h * h);

                C_new[s][i][j] = (C_n[s][i][j] + dt * (-1.5 * H_C + 0.5 * H_C_old
                                                       + Df[s] * lapC
                                                       + ramp * I_S[i][j] * Cs[s][i][j] / dtau)) /
                                 (1. + ramp * I_S[i][j] * dt / dtau);
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