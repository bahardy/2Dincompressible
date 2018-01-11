#ifndef MEMORY_H_INCLUDED
#define MEMORY_H_INCLUDED

#define MEMORY_ALLOCATION \
double** A = make2DDoubleArray(m,n);\
double** Aold = make2DDoubleArray(m,n);\
double** B = make2DDoubleArray(m,n);\
double** Bold = make2DDoubleArray(m,n);\
double*** C = make3DDoubleArray(Ns,m,n);\
double*** Cold = make3DDoubleArray(Ns,m,n); \
double*** Cs = make3DDoubleArray(Ns,m,n); \
double** Cp = make2DDoubleArray(Np,Ns); \
double** C_T = make2DDoubleArray(m,n);\
double** C_Told = make2DDoubleArray(m,n);\
double*** C_C = make3DDoubleArray(Ns,m,n); \
double*** C_Cold = make3DDoubleArray(Ns,m,n); \
double** F = make2DDoubleArray(Np,3);\
double** G = make2DDoubleArray(Np,3);\
double*** Ip_S = make3DDoubleArray(Np,m,n);\
double*** Ip_U = make3DDoubleArray(Np,m,n);\
double*** Ip_V = make3DDoubleArray(Np,m,n);\
double** I_S = make2DDoubleArray(m,n);\
double** coloring = make2DDoubleArray(m,n);\
double** I_U = make2DDoubleArray(m,n);\
double** I_V = make2DDoubleArray(m,n);\
double** M = make2DDoubleArray(Np,3);\
double** P = make2DDoubleArray(m,n);\
double** phi = make2DDoubleArray(m,n);\
double** Qp = make2DDoubleArray(Np,3);\
double*** Qmp = make3DDoubleArray(Np, Ns, 3); \
double** R = make2DDoubleArray(m,n);\
double** T = make2DDoubleArray(m,n);\
double** Told = make2DDoubleArray(m,n);\
double* Tp = make1DDoubleArray(Np);\
double** Ts = make2DDoubleArray(m,n); \
double** U = make2DDoubleArray(m,n);\
double** Up = make2DDoubleArray(Np,4);\
double** Us = make2DDoubleArray(m,n);\
double** Ustar = make2DDoubleArray(m,n);\
double** V = make2DDoubleArray(m,n);\
double** Vp = make2DDoubleArray(Np,4);\
double** Vs = make2DDoubleArray(m,n);\
double** Vstar = make2DDoubleArray(m,n);\
double** Wp = make2DDoubleArray(Np,4);\

#define FREE_MEMORY \
for(int s=0; s<Ns; s++){\
    for(int i=0; i<m; i++){ \
        free(C[s][i]); free(Cold[s][i]); free(Cs[s][i]); \
        free(C_C[s][i]); free(C_Cold[s][i]);\
    }\
}\
for(int s=0; s<Ns; s++){\
    free(C[s]); free(Cold[s]); free(Cs[s]);\
    free(C_C[s]); free(C_Cold[s]);\
}\
for(int k=0; k<Np; k++){ \
    for(int i=0; i<m; i++){ \
        free(Ip_S[k][i]); free(Ip_U[k][i]); free(Ip_V[k][i]);\
    }\
    for(int s=0; s<Ns; s++){\
        free(Qmp[k][s]); \
    }\
}\
for(int k=0; k<Np; k++){ \
    free(Ip_S[k]); free(Ip_U[k]); free(Ip_V[k]);\
    free(Up[k]); free(Vp[k]); free(Wp[k]);\
    free(F[k]); free(G[k]); free(M[k]); \
    free(Qp[k]); free(Qmp[k]);\
    free(Cp[k]);\
}\
\
for(int i=0; i<m; i++){\
    free(I_S[i]); free(I_U[i]); free(I_V[i]); free(coloring[i]);  \
    free(U[i]); free(Ustar[i]); free(V[i]); free(Vstar[i]);\
    free(A[i]); free(Aold[i]); free(B[i]); free(Bold[i]);\
    free(C_T[i]); free(C_Told[i]); free(Told[i]); \
    free(Us[i]); free(Vs[i]); free(Ts[i]);\
    free(P[i]); free(phi[i]); free(T[i]);\
    free(R[i]);\
}\
\
free(F); free(G); free(M);\
free(Ip_S); free(Ip_U); free(Ip_V);\
free(I_S); free(I_U); free(I_V); free(coloring); \
free(U); free(V); free(Ustar); free(Vstar);\
free(A); free(Aold); free(B); free(Bold);\
free(C); free(Cold); free(Cs); free(Cp); \
free(C_C); free(C_Cold);\
free(C_T); free(C_Told); \
free(Df); \
free(Us); free(Vs); \
free(Up); free(Vp); free(Wp); free(Tp); \
free(P); free(phi); free(T); free(Ts); free(Told);  \
free(Qp); free(Qmp);\
free(R);\
free(xg); free(yg); free(Sp); free(rp); free(dp); free(II);


#endif // MEMORY_H_INCLUDED
