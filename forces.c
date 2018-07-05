//
// Created by Baptiste Hardy on 8/06/18.
//

#include "main.h"
#include "forces.h"

int integrate_penalization(Data *data, double* surf, int k)
{
    // Integral terms
    double** F = data->F;
    double** G = data->G;
    double** M = data->Mz;
    double** QQ = data-> QQ;
    double*** PP = data-> PP;

    double*** Ip_U = data->Ip_U;
    double*** Ip_V = data->Ip_V;
    double*** Ip_S = data->Ip_S;
    double** u_n = data->u_star;
    double** v_n = data->v_star;
    double*** C_n = data->C_n;
    double** T_n = data->T_n;
    double** u_s = data->u_s;
    double** v_s = data->v_s;
    double** Ts = data->Ts;
    double*** Cs = data->Cs;
    double** xg = data->xg;
    double** yg = data->yg;
    double* rp = data->rp;
    int Ns = data->Ns;
    int m = data->m;
    int n = data->n;
    double h = data->h;
    double dtau = data->dtau;

    /* Force along x-direction */
    F[k][0] = F[k][1]; /* n-2*/
    F[k][1] = F[k][2]; /* n-1*/
    //F[k][2] = 0.; /* n*/

    /* Force along y-direction */
    G[k][0] = G[k][1];
    G[k][1] = G[k][2];
    //G[k][2] = 0.;

    /* Moment along z-direction */\
    M[k][0] = M[k][1];
    M[k][1] = M[k][2];
    //M[k][2] = 0.;

#ifdef TEMP
    /* Particle heat balance  */
    QQ[k][0] = QQ[k][1]; /* n-2*/
    QQ[k][1] = QQ[k][2]; /* n-1*/
    //QQ[k][2] = 0.; /* n*/

    for(int s=0; s<Ns; s++){
        PP[k][s][0] = PP[k][s][1];
        PP[k][s][1] = PP[k][s][2];
        //  PP[k][s][2] = 0.;
    }
#endif

//    int startX = (int) floor((xg[k] - rp[k]) / h);
//    PetscPrintf(PETSC_COMM_WORLD, "startX = %d \t", startX);
//    int endX = (int) ceil((xg[k]+rp[k])/h);
//    PetscPrintf(PETSC_COMM_WORLD,"endX = %d \t", endX);
//
//    int startY = (int) floor((yg[k]-rp[k])/h);
//    PetscPrintf(PETSC_COMM_WORLD,"startY = %d \t", startY);
//    int endY = (int) ceil((yg[k]+rp[k])/h);
//    PetscPrintf(PETSC_COMM_WORLD,"endY = %d \t \n", endY);

//    if(startY <= 1||endY >= n-1){
//        PetscPrintf(PETSC_COMM_WORLD,"Wall collision! \n");
//        return 1;
//    }

//    if(endX >= m-1){
//        PetscPrintf(PETSC_COMM_WORLD,"Particle leaves the channel! \n");
//        return 1;
//    }

    double Fint, Gint, Mint, Qint, *Qmint, sint;
    Fint = 0.;
    Gint = 0.;
    Mint = 0.;
    Qint = 0.;
    Qmint = make1DDoubleArray(Ns);
    sint = 0.;

    double h2;
    h2 = h*h;
    double yU, xV, f, g, q, qm;


    //for(int i=startX; i<=endX; i++){//
    for(int i = 0; i<m; i++){
        xV = (i-0.5)*h;
        //for(int j=startY; j<=endY; j++){
        for(int j=1; j<n-1; j++) {
            yU = (j-0.5)*h;
            f = -Ip_U[k][i][j]*(u_n[i][j]-u_s[i][j]);
            g = -Ip_V[k][i][j]*(v_n[i][j]-v_s[i][j]);
#ifdef TEMP
            q = -Ip_S[k][i][j]*(T_n[i][j]-Ts[i][j]);
            Qint += q;
            for(int s=0; s<Ns; s++){
                qm = -Ip_S[k][i][j]*(C_n[s][i][j]-Cs[s][i][j]);
                Qmint[s] += qm;
            }
#endif
            sint += Ip_S[k][i][j]*h*h;

            Fint += f; /* units : m/s */
            Gint += g; /* units : m/s */
            Mint += ((xV-xg[k][1])*g-(yU-yg[k][1])*f);/* units: m^2/s */
        }
    }
    Fint *= h2/dtau; /* units : m^3/s; */
    Gint *= h2/dtau;
    Mint *= h2/dtau;
    Qint *= h2/dtau; /* units : K*m^2/s */

    F[k][2] = -Fint;
    G[k][2] = -Gint;
    M[k][2] = -Mint;
    QQ[k][2] = -Qint;
    *surf = sint;

#ifdef TEMP
    for(int s=0; s<Ns; s++){
        PP[k][s][2] = -Qmint[s]*h2/dtau;
    }
#endif
    free(Qmint);
    PetscPrintf(PETSC_COMM_WORLD, "Particle surface is %f\n", *surf);
    return 0;
}

void integrate_penalization_periodic(Data *data, double* surf, int k)
{
    // Integral terms
    double** F = data->F;
    double** G = data->G;
    double** M = data->Mz;
    double** QQ = data-> QQ;
    double*** PP = data-> PP;

    double*** Ip_U = data->Ip_U;
    double*** Ip_V = data->Ip_V;
    double*** Ip_S = data->Ip_S;
    double** u_n = data->u_star;
    double** v_n = data->v_star;
    double*** C_n = data->C_n;
    double** T_n = data->T_n;
    double** u_s = data->u_s;
    double** v_s = data->v_s;
    double** Ts = data->Ts;
    double*** Cs = data->Cs;
    double* d = make1DDoubleArray(3);
    double dx_min;

    int Ns = data->Ns;
    int m = data->m;
    int n = data->n;
    double L = data->L;
    double h = data->h;
    double dtau = data->dtau;


#ifdef TEMP
    /* Particle heat balance  */
    QQ[k][0] = QQ[k][1]; /* n-2*/
    QQ[k][1] = QQ[k][2]; /* n-1*/
    //QQ[k][2] = 0.; /* n*/

    for(int s=0; s<Ns; s++){
        PP[k][s][0] = PP[k][s][1];
        PP[k][s][1] = PP[k][s][2];
        //  PP[k][s][2] = 0.;
    }
#endif

    double Fint, Gint, Mint, Qint, *Qmint, sint;
    Fint = 0.;
    Gint = 0.;
    Mint = 0.;
    Qint = 0.;
    Qmint = make1DDoubleArray(Ns);
    sint = 0.;

    double h2 = h*h;
    double yU, xV, f, g, q, qm;

    double xG = fmod(data->xg[k][1], L);
    double YG = data->yg[k][1];

    int i_min;

    for(int i = 0; i<m; i++){
        xV = (i+0.5)*h;
        d[0] = xV - (xG-L);
        d[1] = xV - xG;
        d[2] = xV - (xG+L);
        i_min = min_abs(d,3);
        dx_min = d[i_min];

        for(int j=1; j<n-1; j++) {
            yU = (j-0.5)*h;
            f = -Ip_U[k][i][j]*(u_n[i][j]-u_s[i][j]);
            g = -Ip_V[k][i][j]*(v_n[i][j]-v_s[i][j]);
#ifdef TEMP
            q = -Ip_S[k][i][j]*(T_n[i][j]-Ts[i][j]);
            Qint += q;
            for(int s=0; s<Ns; s++){
                qm = -Ip_S[k][i][j]*(C_n[s][i][j]-Cs[s][i][j]);
                Qmint[s] += qm;
            }
#endif
            sint += Ip_S[k][i][j]*h*h;

            Fint += f; /* units : m/s */
            Gint += g; /* units : m/s */
            Mint += dx_min*g-(yU-YG)*f;/* units: m^2/s */
        }
    }
    Fint *= h2/dtau; /* units : m^3/s; */
    Gint *= h2/dtau;
    Mint *= h2/dtau;
    Qint *= h2/dtau; /* units : K*m^2/s */

    F[k][2] = -Fint;
    G[k][2] = -Gint;
    M[k][2] = -Mint;
    QQ[k][2] = -Qint;
    *surf = sint;

#ifdef TEMP
    for(int s=0; s<Ns; s++){
        PP[k][s][2] = -Qmint[s]*h2/dtau;
    }
#endif
    free(Qmint);
    free(d);
}

void compute_forces_fluxes(Data* data, int k)
{

    double dt = data->dt;

    double* Fx = data->Fx;
    double* Fy = data->Fy;
    double* Tz = data->Tz;
    double* Q = data->Q;
    double** Qm = data->Qm;

    double** F = data->F;
    double Fbis = 0;
    double Fter = 0;
    double** G = data->G;
    double** M = data->Mz;
    double** QQ = data->QQ;
    double*** PP = data->PP;

    double** Up = data->Up;
    double** Vp = data->Vp;
    double** Omega_p = data->Omega_p;
    double* dudt = data->dudt;
    double* dvdt = data->dvdt;
    double* domegadt = data->domegadt;
    double* dTdt = data->dTdt;

    double* Sp = data->Sp;
    double* J = data->J;

    double rho_f = data->rho_f;
    double rho_p = data->rho_p;
    double rho_r = data->rho_r;

    double cf = data->cf;

    PetscPrintf(PETSC_COMM_WORLD,"F integration = %1.6e \n", k+1, F[k][2]);
    PetscPrintf(PETSC_COMM_WORLD,"G integration = %1.6e \n", k+1, G[k][2]);
    PetscPrintf(PETSC_COMM_WORLD,"M integration = %1.6e \n", k+1, M[k][2]);

    dudt[k] = (Up[k][2]-Up[k][0])/(2*dt);
    dvdt[k] = (Vp[k][2]-Vp[k][0])/(2*dt);
    domegadt[k] = (Omega_p[k][2]-Omega_p[k][0])/(2*dt);

    PetscPrintf(PETSC_COMM_WORLD,"dudt = %1.6e \n", k+1, dudt[k]);
    PetscPrintf(PETSC_COMM_WORLD,"dvdt = %1.6e \n", k+1, dvdt[k]);
    PetscPrintf(PETSC_COMM_WORLD,"domegadt = %1.6e \n", k+1, domegadt[k]);

    Fx[k] = rho_f*(Sp[k]*dudt[k] + F[k][2]);
    Fbis = rho_p*Sp[k]*dudt[k];
    Fter = rho_f*(rho_r/(rho_r -1))*F[k][2];
    Fy[k] = rho_f*(Sp[k]*dvdt[k] + G[k][2]);
    Tz[k] = rho_f*(J[k]*domegadt[k] + M[k][2]);
    Q[k] = rho_f*cf*(Sp[k]*dTdt[k] + QQ[k][2]);
    Qm[k][0] = PP[k][0][2];


    PetscPrintf(PETSC_COMM_WORLD,"Hydrodynamic force along -x dir on particle %d = %1.6e [N/m] OR = %1.6e [N/m] OR = %1.6e [N/m] \n", k+1, Fx[k], Fbis, Fter);
    PetscPrintf(PETSC_COMM_WORLD,"Hydrodynamic force along -y dir on particle %d = %1.6e [N/m] \n", k+1, Fy[k]);
    PetscPrintf(PETSC_COMM_WORLD,"Torque on particle %d = %1.6e [N]  \n", k+1, Tz[k]);
    //PetscPrintf(PETSC_COMM_WORLD,"Heat flux on particle %d = %1.6e [W/m] \n", k+1, Q[k]);
    //PetscPrintf(PETSC_COMM_WORLD,"Molar flux of A on particle %d = %1.6e [mol/(m.s)] \n", k+1, Qm[k][0]);
}

void get_tau(Data* data)
{
    int m = data->m;
    int n = data->n;
    int i, j;
    double nu = data->nu;
    double dx = data->h;
    double dy = data->h;
    double** u_n = data->u_n;
    double** v_n = data->v_n;
    double** tau_xx = data->tau_xx;
    double** tau_yy = data->tau_yy;
    double** tau_xy = data->tau_xy;

    double dudx, dudy, dvdy, dvdx;

    /** Diagonal components **/
    for(i = 1; i<m-1; i++) {
        for(j=1; j<n-1; j++){
            dudx = (u_n[i][j] - u_n[i-1][j])/dx;
            dvdy = (v_n[i][j] - v_n[i][j-1])/dy;
            tau_xx[i][j] = 2*nu*dudx;
            tau_yy[i][j] = 2*nu*dvdy;
        }
    }
    /** Off-diagonal components **/
    for(i=0; i<m-1; i++){
        for(j=0; i<n-1; i++) {
            dudy = (u_n[i][j+1] - u_n[i][j])/dy;
            dvdx = (v_n[i+1][j] - v_n[i][j])/dx;
            tau_xy[i][j] = nu*(dudy + dvdx);
        }
    }
}

void get_tau_periodic(Data* data)
{
    int m = data->m;
    int n = data->n;
    int i, j;
    double nu = data->nu;
    double dx = data->h;
    double dy = data->h;
    double** u_n = data->u_n;
    double** v_n = data->v_n;
    double** tau_xx = data->tau_xx;
    double** tau_yy = data->tau_yy;
    double** tau_xy = data->tau_xy;

    double dudx, dudy, dvdy, dvdx;

    /** Diagonal components **/
    for(i = 0; i<m; i++) {
        for(j=1; j<n-1; j++){
            dudx = (u_n[i][j] - u_n[(i-1+m)%m][j])/dx;
            dvdy = (v_n[i][j] - v_n[i][j-1])/dy;
            tau_xx[i][j] = 2*nu*dudx;
            tau_yy[i][j] = 2*nu*dvdy;
        }
    }
    /** Off-diagonal components **/
    for(i=0; i<m; i++){
        for(j=0; i<n-1; i++) {
            dudy = (u_n[i][j+1] - u_n[i][j])/dy;
            dvdx = (v_n[(i+1+m)%m][j] - v_n[i][j])/dx;
            tau_xy[i][j] = nu*(dudy + dvdx);
        }
    }
}

void compute_forces_NOCA(Data* data, FILE* file, int I1, int I2, int J1, int J2)
{
    int i, j;
    double ILx, ITx, IRx, IBx;
    double ILy, ITy, IRy, IBy;
    ILx = 0, ITx = 0, IRx = 0, IBx = 0;
    ILy = 0, ITy = 0, IRy = 0, IBy = 0;

    double IU, IV;
    IU = 0, IV = 0;

    double h = data->h;
    double** u_n = data->u_n;
    double** v_n = data->v_n;
    double** chi_U = data->I_U;
    double** chi_V = data->I_V;
    double** P = data->P;
    double** tau_xx = data->tau_xx;
    double** tau_xy = data->tau_xy;
    double** tau_yy = data->tau_xy;
    double I_cs_u = 0;
    double I_cs_v = 0;
    double u_av, v_av;

    for(j = J1; j<J2; j++)
    {
        ILx += (u_n[I1][j]*u_n[I1][j] + .5*(P[I1][j] + P[I1+1][j]) - .5*(tau_xx[I1][j] + tau_xx[I1+1][j]))*h;
        v_av = .25*(v_n[I1][j] + v_n[I1+1][j] + v_n[I1][j-1] + v_n[I1+1][j-1]);
        ILy += (u_n[I1][j]*v_av - .5*(tau_xy[I1][j] + tau_xy[I1][j-1]))*h;

        IRx += (-u_n[I2][j]*u_n[I2][j] - .5*(P[I2][j] + P[I2+1][j]) + .5*(tau_xx[I2][j] + tau_xx[I2+1][j]))*h;
        v_av = .25*(v_n[I2][j] + v_n[I2+1][j] + v_n[I2][j-1] + v_n[I2+1][j-1]);
        IRy += (-u_n[I2][j]*v_av + .5*(tau_xy[I2][j] + tau_xy[I2][j-1]))*h;

    }

    for(i = I1; i<I2; i++)
    {
        u_av = .25*(u_n[i][J2] + u_n[i][J2+1] + u_n[i-1][J2] + u_n[i-1][J2+1]);
        ITx += (-v_n[i][J2]*u_av + .5*(tau_xy[i][J2] + tau_xy[i-1][J2]) )*h;
        ITy += (-v_n[i][J2]*v_n[i][J2] - .5*(P[i][J2]+P[i][J2+1]) + .5*(tau_yy[i][J2] + tau_yy[i][J2+1]))*h;

        u_av = .25*(u_n[i][J1] + u_n[i][J1+1] + u_n[i-1][J1] + u_n[i-1][J1+1]);
        IBx += (v_n[i][J1]*u_av - .5*(tau_xy[i][J1] + tau_xy[i-1][J1]) )*h;
        IBy += (v_n[i][J1]*v_n[i][J1] + .5*(P[i][J1]+P[i][J1+1]) - .5*(tau_yy[i][J1] + tau_yy[i][J1+1]))*h;
    }

    I_cs_u = ILx + ITx + IRx + IBx;
    I_cs_v = ILy + ITy + IRy + IBy;


    for(i=I1; i<I2; i++) {
        for (j = J1 + 1; j < J2; j++) {
            IU += (1-chi_U[i][j]) * u_n[i][j] * h * h;
        }
    }

    for (i = I1+1; i<I2; i++) {
        for(j = J1; j <J2; j++) {
            IV += (1-chi_V[i][j])* v_n[i][j] * h * h;
        }
    }

    fprintf(file, "%3.10e \t %3.10e \t %3.10e \t %3.10e \n", IU, I_cs_u, IV, I_cs_v);
    fflush(file);

}

int min_abs(double* array, int arraySize)
{
    double min = INFINITY;
    int index = -1;
    for(int i = 0; i < arraySize; i++) {
        if (fabs(array[i]) < min) {
            min = fabs(array[i]);
            index = i;
        }
    }
    return index;
}