//
// Created by Baptiste Hardy on 4/05/18.
//

#include "main.h"
#include "fields_memory.h"

void allocate_fields(Data* data)
{
    int m = data->m;
    int n = data->n;
    int Np = data->Np;
    int Ns = data->Ns;

    data->coloring = make2DDoubleArray(m,n);
    data->CFL_array = make2DDoubleArray(m,n);
    data->C_n = make3DDoubleArray(Ns,m,n);
    data->C0 = make1DDoubleArray(Ns); // inlet concentration
    data->C_n_1 = make3DDoubleArray(Ns,m,n);
    data->Cs = make3DDoubleArray(Ns,m,n);
    data->Cp = make2DDoubleArray(Np,Ns);
    data->dudt = make1DDoubleArray(Np);
    data->dvdt = make1DDoubleArray(Np);
    data->domegadt = make1DDoubleArray(Np);
    data->dTdt = make1DDoubleArray(Np);
    data->dCdt = make2DDoubleArray(Np,Ns);
    data->dp = make1DDoubleArray(Np);

    data->Fx = make1DDoubleArray(Np);
    data->Fy = make1DDoubleArray(Np);
    data->Tz = make1DDoubleArray(Np);
    data->Fx_coll = make2DDoubleArray(Np,3);
    data->Fy_coll = make2DDoubleArray(Np,3);

    data->F = make2DDoubleArray(Np,3);
    data->G = make2DDoubleArray(Np,3);
    data->Mz = make2DDoubleArray(Np,3);

    data->chi_S = make3DDoubleArray(Np,m,n);
    data->chi_U = make3DDoubleArray(Np,m,n);
    data->chi_V = make3DDoubleArray(Np,m,n);
    data->Ip_S = make3DDoubleArray(Np,m,n);
    data->Ip_U = make3DDoubleArray(Np,m,n);
    data->Ip_V = make3DDoubleArray(Np,m,n);
    data->I_S = make2DDoubleArray(m,n);
    data->I_U = make2DDoubleArray(m,n);
    data->I_V = make2DDoubleArray(m,n);

    data->J = make1DDoubleArray(Np);
    data->phi = make2DDoubleArray(m,n);
    data->Qm = make2DDoubleArray(Np,Ns);
    data->P = make2DDoubleArray(m,n);
    data->PP = make3DDoubleArray(Np, Ns, 3);
    data->Q = make1DDoubleArray(Np);
    data->QQ = make2DDoubleArray(Np,3);
    data->Qr = make2DDoubleArray(Np,3);
    data->rp = make1DDoubleArray(Np);
    data->Reh = make2DDoubleArray(m,n);
    data->Reh_omega = make2DDoubleArray(m,n);
    data->Sp = make1DDoubleArray(Np);

    data->H_u_n_1 = make2DDoubleArray(m,n);
    data->H_v_n_1 = make2DDoubleArray(m,n);
    data->H_T_n_1 = make2DDoubleArray(m,n);
    data->H_C_n_1 = make3DDoubleArray(Ns, m, n);

    data->T_n = make2DDoubleArray(m,n);
    data->T_n_1 = make2DDoubleArray(m,n);
    data->Tp = make1DDoubleArray(Np);
    data->Ts = make2DDoubleArray(m,n);

    data->u_n = make2DDoubleArray(m,n);
    data->u_n_1 = make2DDoubleArray(m,n);
    data->u_s = make2DDoubleArray(m,n);
    data->u_star = make2DDoubleArray(m,n);

    data->omega = make2DDoubleArray(m,n);

    data->v_n = make2DDoubleArray(m,n);
    data->v_n_1 = make2DDoubleArray(m,n);
    data->v_s = make2DDoubleArray(m,n);
    data->v_star = make2DDoubleArray(m,n);

    data->tau_xx = make2DDoubleArray(m,n);
    data->tau_yy = make2DDoubleArray(m,n);
    data->tau_xy = make2DDoubleArray(m,n);

    data->xg = make2DDoubleArray(Np,3);
    data->yg = make2DDoubleArray(Np,3);
    data->theta = make2DDoubleArray(Np,3);

    data->Up = make2DDoubleArray(Np,3);
    data->Vp = make2DDoubleArray(Np,3);
    data->Omega_p = make2DDoubleArray(Np,3);


}

void free_fields(Data* data)
{
    int m = data->m;
    int Np = data->Np;
    int Ns = data->Ns;

    /* Free memory */
    free2Darray(data->xg,Np), free2Darray(data->yg,Np), free2Darray(data->theta,Np);
    free(data->dp), free(data->rp), free(data->Sp), free(data->J);
    free(data->dudt), free(data->dvdt), free(data->domegadt), free(data->dTdt); free2Darray(data->dCdt, Np);
    free(data->Fx), free(data->Fy), free(data->Tz), free(data->Q), free2Darray(data->Qm, Np);
    free2Darray(data->Fx_coll, Np), free2Darray(data->Fy_coll, Np);
    free2Darray(data->u_n,m), free2Darray(data->u_n_1,m), free2Darray(data->u_star, m), free2Darray(data->u_s, m);
    free2Darray(data->H_u_n_1, m) , free2Darray(data->H_v_n_1, m); //free2Darray(data->H_T_n_1, m), free3Darray(data->H_Y_n_1, Ns, m);
    free2Darray(data->v_n,m), free2Darray(data->v_n_1,m), free2Darray(data->v_star,m), free2Darray(data->v_s,m);
    free2Darray(data->omega, m); free2Darray(data->Reh,m); free2Darray(data->Reh_omega,m); free2Darray(data->CFL_array, m);
    free2Darray(data->P,m), free2Darray(data->phi, m);
    free2Darray(data->T_n,m),  free2Darray(data->T_n_1,m), free2Darray(data->Ts, m);
    free3Darray(data->C_n, Ns, m), free3Darray(data->C_n_1, Ns, m), free3Darray(data->Cs, Ns, m), free(data->C0);
    free2Darray(data->Up,Np), free2Darray(data->Vp, Np), free2Darray(data->Omega_p,Np), free(data->Tp), free2Darray(data->Cp, Np);
    free2Darray(data->F, Np), free2Darray(data->G, Np), free2Darray(data->Mz, Np), free2Darray(data->QQ, Np), free3Darray(data->PP, Np,Ns), free2Darray(data->Qr, Np);
    free2Darray(data->I_S, m), free2Darray(data->I_U, m), free2Darray(data->I_V, m), free2Darray(data->coloring, m);
    free3Darray(data->Ip_S, Np,m), free3Darray(data->Ip_U, Np,m), free3Darray(data->Ip_V, Np, m);
    free3Darray(data->chi_S, Np, m), free3Darray(data->chi_U, Np, m), free3Darray(data->chi_V, Np, m);
    free2Darray(data->tau_xx, m), free2Darray(data->tau_yy, m),free2Darray(data->tau_xy, m);
}

