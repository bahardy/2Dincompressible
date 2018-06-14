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
    double* xg = data->xg;
    double* yg = data->yg;
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

    int startX = (int) floor((xg[k] - rp[k]) / h);
    PetscPrintf(PETSC_COMM_WORLD, "startX = %d \t", startX);
    int endX = (int) ceil((xg[k]+rp[k])/h);
    PetscPrintf(PETSC_COMM_WORLD,"endX = %d \t", endX);

    int startY = (int) floor((yg[k]-rp[k])/h);
    PetscPrintf(PETSC_COMM_WORLD,"startY = %d \t", startY);
    int endY = (int) ceil((yg[k]+rp[k])/h);
    PetscPrintf(PETSC_COMM_WORLD,"endY = %d \t \n", endY);

    if(startY <= 1||endY >= n-1){
        PetscPrintf(PETSC_COMM_WORLD,"Wall collision! \n");
        return 1;
    }

    if(endX >= m-1){
        PetscPrintf(PETSC_COMM_WORLD,"Particle leaves the channel! \n");
        return 1;
    }

    double Fint, Gint, Mint, Qint, *Qmint, sint;
    Fint = 0.;
    Gint = 0.;
    Mint = 0.;
    Qint = 0.;
    Qmint = make1DDoubleArray(Ns);
    sint = 0.;

    double h2;
    h2 = h*h;
    double yU, xV, f, g, q, *qm;
    qm = make1DDoubleArray(Ns);


    //for(int i=startX; i<=endX; i++){//
    for(int i = 0; i<m; i++){
        xV = (i-0.5)*h;
        //for(int j=startY; j<=endY; j++){
        for(int j=0; j<n; j++) {
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
            Mint += ((xV-xg[k])*g-(yU-yg[k])*f);/* units: m^2/s */
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
    free(qm);
    PetscPrintf(PETSC_COMM_WORLD, "Particle surface is %f\n", *surf);
    return 0;
}

void compute_forces_fluxes(Data* data, int k)
{

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

    PetscPrintf(PETSC_COMM_WORLD,"\n F integration = %1.6e \n", k+1, F[k][2]);
    PetscPrintf(PETSC_COMM_WORLD,"G integration = %1.6e \n", k+1, G[k][2]);
    PetscPrintf(PETSC_COMM_WORLD,"M integration = %1.6e \n", k+1, M[k][2]);

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
    PetscPrintf(PETSC_COMM_WORLD,"Heat flux on particle %d = %1.6e [W/m] \n", k+1, Q[k]);
    PetscPrintf(PETSC_COMM_WORLD,"Molar flux of A on particle %d = %1.6e [mol/(m.s)] \n", k+1, Qm[k][0]);
}
