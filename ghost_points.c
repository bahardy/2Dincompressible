//
// Created by Baptiste Hardy on 13/08/18.
//

#include "ghost_points.h"

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
        /** IMPOSED TEMP **/
        double T_left = 1;
        double T_right = 0;
        T_n[i][0] = -0.2*(T_n[i][3]-5.*T_n[i][2]+15.*T_n[i][1]-16.*T_left);
        T_n[i][n-1] = -0.2*(T_n[i][n-4]-5.*T_n[i][n-3]+15.*T_n[i][n-2]-16.*T_right);

        /* Walls : adiabatic: dTdn = 0, no mass flux */
        //T_n[i][0] = T_n[i][1];
        //T_n[i][n-1] = T_n[i][n-2];
        for (int s=0; s<Ns; s++){
//            C[s][i][0] = -0.2*(C[s][i][3]-5.*C[s][i][2]+15.*C[s][i][1]-16.*C0[0]);
//            C[s][i][n-1] = -0.2*(C[s][i][n-4]-5.*C[s][i][n-3]+15.*C[s][i][n-2]-16.*C0[0]);

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
        //v_n[m-1][j] = v_n[m-2][j];
        /* Closed cavity */
        v_n[m-1][j] = -0.2*(v_n[m-4][j] - 5.*v_n[m-3][j] + 15.*v_n[m-2][j]);

#ifdef TEMP
        /* On T_n and C */
        /* Inflow : T_n uniform  */
        /* Outflow : We cancel axial dispersion d2T/dx2 = 0; d2C/dx2 = 0; */

        /** ADIABATIC **/
        T_n[0][j] = T_n[1][j];
        T_n[m-1][j] = T_n[m-2][j];

        //T_n[0][j] = -0.2*(T_n[3][j]-5.*T_n[2][j]+15.*T_n[1][j]-16.*T0);
        //T_n[m-1][j] = (7.*T_n[m-2][j]-5.*T_n[m-3][j]+T_n[m-4][j])/3.;

        for(int s=0; s<Ns; s++){
            C[s][0][j] = -0.2*(C[s][3][j]-5.*C[s][2][j]+15.*C[s][1][j]-16.*C0[s]);
            C[s][m-1][j] = (7.*C[s][m-2][j]-5.*C[s][m-3][j]+C[s][m-4][j])/3.;
            //C[s][m-1][j] = -0.2*(C[s][m-4][j]-5.*C[s][m-3][j]+15.*C[s][m-2][j]-16.*C0[s]);
        }

#endif

    }
}
