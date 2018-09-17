//
// Created by Baptiste Hardy on 18/06/18.
//

#include "main.h"
#include "collision.h"

int sign(double x)
{
    int sign = (x > 0) - (x < 0);
    return sign;
}

void collision(Data* data)
{

    int k1, k2;
    double x1, y1, R1, x2, y2, R2;
    double z = data->h;
    double d12, d_wall;
    double d = data->d;
    double** xg = data->xg;
    double** yg = data->yg;
    double* rp = data->rp;
    double** Up = data->Up;
    double** Vp = data->Vp;
    double** Omega_p = data->Omega_p;
    double** Fx_coll = data->Fx_coll;
    double** Fy_coll = data->Fy_coll;
    double** M_coll = data->M_coll;

    int Np = data->Np;
    double rho_s = data->rho_s;
    double rho_f = data->rho_f;
    double g = data->g;
    double S = data->Sp[0];
    double L = data->L;


    double c11 = S*(rho_s - rho_f)*g;
    double c12 = c11;
    double ep = data->ep;
    double Ep = data->Ep;
    double ew = data->ew;//1e-6;
    double Ew = data->Ew; //1e-8;

    double uct, vct, vn, nx, ny, vct_norm;
    double u1, u2, v1, v2, W1, W2;
    double F_coll_n = 0;
    double F_coll_t;
    double Fnx, Fny, Ftx, Fty;
    double tx, ty;
    double mu_c = 0.8;

    for(k1=0; k1<Np; k1++) {
        x1 = fmod(xg[k1][2],L);
        y1 = yg[k1][2];
        R1 = rp[k1];
        u1 = Up[k1][2];
        v1 = Vp[k1][2];
        W1 = Omega_p[k1][2];
        for(k2=k1+1; k2<Np; k2++) {
            x2 = fmod(xg[k2][2],L);
            y2 = yg[k2][2];
            R2 = rp[k2];
            d12 = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
            d12 = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
            nx = -(x1 - x2)/d12;
            ny = -(y1 - y2)/d12;
            Fnx = 0;
            Fny = 0;
            Ftx = 0;
            Fty = 0;

            if (d12 <=(R1 + R2 + z))
            {
                // AS SOON AS THERE IS A COLLISION
                vn = (u1 - u2)*nx + (v1 - v2)*ny; // relative normal velocity
                uct = (u1-u2) - vn*nx - (R1*W1 + R2*W2)*ny; // x-compo of the tangential relative velocity
                vct = (v1-v2) - vn*ny + (R1*W1 + R2*W2)*nx; // y-compo of the tangential relative velocity

                vct_norm = sqrt(uct*uct + vct*vct);
                tx = uct/vct_norm;
                ty = vct/vct_norm;

                if (d12 > (R1 + R2)){
                    F_coll_n = (c12/ep)*pow((d12 - (R1 + R2 + z))/z , 2.);
                }
                else {
                    F_coll_n = ((c12/ep)*pow((d12 - (R1 + R2 + z))/z , 2.) + (c12/Ep)*(R1+R2-d12)/z);
                }
                Fnx = -F_coll_n*nx;
                Fny = -F_coll_n*ny;

                F_coll_t = mu_c*F_coll_n; // Ft = -mu*Fn*t
                Ftx = -F_coll_t*tx;
                Fty = -F_coll_t*ty;


            }

            Fx_coll[k1][2] += Fnx + Ftx;
            Fx_coll[k2][2] += -Fnx - Ftx;
            Fy_coll[k1][2] += Fny + Fty;
            Fy_coll[k2][2] += -Fny - Fty;

            M_coll[k1][2] += R1*(nx*Fty - ny*Ftx);
            M_coll[k2][2] += R2*(nx*Fty - ny*Ftx);

        }
        /** Check for wall collisions **/

        d_wall = fmin(y1, d-y1);

        if(2*d_wall < 2*R1)
        {
            Fy_coll[k1][2] += ( (c11/ew)*pow((2*d_wall - (2*R1 + z))/z, 2.) + (c12/Ew)*(2*R1 - 2*d_wall)/z )*sign(d/2. - y1);
        }
        else if (2*d_wall < 2*R1 + z)
        {
            Fy_coll[k1][2] += (c11/ew)*pow((2*d_wall - (2*R1 + z))/z, 2.)*sign(d/2. - y1);
        }

        d_wall = fmin(x1, L-x1);

        if(2*d_wall < 2*R1)
        {
            Fx_coll[k1][2] += ( (c11/ew)*pow((2*d_wall - (2*R1 + z))/z, 2.) + (c12/Ew)*(2*R1 - 2*d_wall)/z )*sign(L/2. - x1);
        }
        else if (2*d_wall < 2*R1 + z)
        {
            Fx_coll[k1][2] += (c11/ew)*pow((2*d_wall - (2*R1 + z))/z, 2.)*sign(L/2. - x1);
        }

    }
}
