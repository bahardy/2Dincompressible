//
// Created by Baptiste Hardy on 8/03/18.
//

#include <stdio.h>
#include <string.h>
#include "write.h"

void writeData(FILE* fichier_data, Data data)
{
    if (fichier_data != NULL){
        fprintf(fichier_data,"dp \t %f\n",data.Dp);
        fprintf(fichier_data,"ratio d/dp \t %f\n",data.d);
        fprintf(fichier_data,"ratio dp/h \t %f\n",data.Dp/data.h);
        fprintf(fichier_data, "ratio L/d \t %f\n", data.L/data.d);
        fprintf(fichier_data,"d\t %f\n",data.d);
        fprintf(fichier_data,"H\t %f\n",data.H);
        fprintf(fichier_data,"rho_p\t %f\n",data.rho_p);
        fprintf(fichier_data,"rho_f\t %f\n",data.rho_f);
        fprintf(fichier_data,"rho_r\t %f\n",data.rho_r);
        fprintf(fichier_data,"cp_s\t %f\n",data.cp_s);
        fprintf(fichier_data,"cp_f\t %f\n",data.cp_f);
        fprintf(fichier_data,"cr\t %f\n",data.cr);
        fprintf(fichier_data,"Rep\t %f\n",data.Rep);
        fprintf(fichier_data,"Pr\t %f\n",data.Pr);
        fprintf(fichier_data,"Sc\t %f\n",data.Sc);
        fprintf(fichier_data,"Le\t %f\n",data.Le);
        fprintf(fichier_data,"Fr\t %f\n",data.Fr);
        fprintf(fichier_data,"nu\t %1.8f\n",data.nu);
        fprintf(fichier_data,"alpha_SOR\t %1.8f\n",data.alpha_f);
        fprintf(fichier_data,"Ns \t %d\n", data.Ns);
        fprintf(fichier_data,"DfA\t %1.8f\n",data.Df[0]);
        fprintf(fichier_data,"DfB\t %1.8f\n",data.Df[1]);
        fprintf(fichier_data,"CA0\t %f\n",data.CA0);
        fprintf(fichier_data,"CB0\t %f\n",data.CB0);
        fprintf(fichier_data,"dH \t %1.8f\n",data.dH);
        fprintf(fichier_data,"Um\t %f\n", data.u_m);
        fprintf(fichier_data,"Umax\t %f\n", data.u_max);
        fprintf(fichier_data,"dt\t %1.8f\n", data.dt);
        fprintf(fichier_data,"dtau\t %1.8f\n", data.dtau);
        fprintf(fichier_data,"ratio dtau/dt \t %1.8f\n", data.ratio_dtau_dt);
        fprintf(fichier_data,"CFL\t %1.8f\n", data.CFL);
        fprintf(fichier_data,"Fourier \t %1.8f\n", data.r);
        fprintf(fichier_data,"n\t %d\n",data.n);
        fprintf(fichier_data,"m\t %d\n",data.m);
        fprintf(fichier_data,"N\t %d\n",data.N);
        fprintf(fichier_data,"M\t %d\n",data.M);
        fprintf(fichier_data,"h\t %f\n",data.h);
        fprintf(fichier_data,"Kmax\t %d\n", data.Kmax);
        fprintf(fichier_data,"nKmax\t %d\n", data.nKmax);
        fprintf(fichier_data,"Nwrite\t %d\n", data.N_write);
        fprintf(fichier_data,"Twrite\t %d\n", data.T_write);
        fprintf(fichier_data,"Tf\t %f\n", data.Tf);
        fprintf(fichier_data, "Np \t %d \n", data.Np);
        fprintf(fichier_data, "Ns \t %d \n", data.Ns);
        fprintf(fichier_data, "Tm0 \t %f \n", data.T0);
        fprintf(fichier_data, "Tp0 \t %f \n", data.Tp0);
        fprintf(fichier_data, "epsilon/h \t %f \n", data.eps/data.h);
        fprintf(fichier_data, "ep \t %f \n", data.ep);
        fprintf(fichier_data, "Ep \t %f \n", data.Ep);
        fprintf(fichier_data, "ew \t %f \n", data.ew);
        fprintf(fichier_data, "Ew \t %f \n", data.Ew);
        fflush(fichier_data);
    }
    else{\
        printf("An error occured in writeFile : invalid file pointer.\n");\
    }
}

void writeFields(Data* data, char* folder, int it) {
    int m = data->m;
    int n = data->n;

    FILE *fichier_U = NULL;
    FILE *fichier_V = NULL;
    //FILE* fichier_P = NULL;
    FILE *fichier_T = NULL;
    FILE *fichier_YA = NULL;
    FILE *fichier_YB = NULL;
    FILE *fichier_Mask = NULL;
    FILE *fichier_nx = NULL;
    FILE *fichier_ny = NULL;
    FILE *fichier_Vtx = NULL;

    char buffer[10];
    sprintf(buffer, "%d", it);

    char fileMask[30];
    strcpy(fileMask, folder);
    strcat(fileMask, "/Mask");
    strcat(fileMask, "-");
    strcat(fileMask, buffer);
    strcat(fileMask, ".txt");
    fichier_Mask = fopen(fileMask, "w");

    char fileNx[30];
    strcpy(fileNx, folder);
    strcat(fileNx, "/nx");
    strcat(fileNx, "-");
    strcat(fileNx, buffer);
    strcat(fileNx, ".txt");
    fichier_nx = fopen(fileNx, "w");

    char fileNy[30];
    strcpy(fileNy, folder);
    strcat(fileNy, "/ny");
    strcat(fileNy, "-");
    strcat(fileNy, buffer);
    strcat(fileNy, ".txt");
    fichier_ny = fopen(fileNy, "w");

    char fileU[30];
    strcpy(fileU, folder);
    strcat(fileU, "/U");
    strcat(fileU, "-");
    strcat(fileU, buffer);
    strcat(fileU, ".txt");
    fichier_U = fopen(fileU, "w");

    char fileV[30];
    strcpy(fileV, folder);
    strcat(fileV, "/V");
    strcat(fileV, "-");
    strcat(fileV, buffer);
    strcat(fileV, ".txt");
    fichier_V = fopen(fileV, "w");

//    char fileP[30];
//    strcpy(fileP, folder);
//    strcat(fileP, "P");
//    strcat(fileP, "-");
//    strcat(fileP, buffer);
//    strcat(fileP, ".txt");
//    fichier_P = fopen(fileP, "w");

    char fileVtx[30];
    strcpy(fileVtx, folder);
    strcat(fileVtx, "/Vtx");
    strcat(fileVtx, "-");
    strcat(fileVtx, buffer);
    strcat(fileVtx, ".txt");
    fichier_Vtx = fopen(fileVtx, "w");

#ifdef TEMP
    char fileT[30];
    strcpy(fileT, folder);
    strcat(fileT, "/T");
    strcat(fileT, "-");
    strcat(fileT, buffer);
    strcat(fileT, ".txt");
    fichier_T = fopen(fileT, "w");

    char fileYA[30];
    strcpy(fileYA, folder);
    strcat(fileYA, "/CA");
    strcat(fileYA, "-");
    strcat(fileYA, buffer);
    strcat(fileYA, ".txt");
    fichier_YA = fopen(fileYA, "w");

//    char fileYB[30];
//    strcpy(fileYB, folder);
//    strcat(fileYB, "/CB");
//    strcat(fileYB, "-");
//    strcat(fileYB, buffer);
//    strcat(fileYB, ".txt");
//    fichier_YB = fopen(fileYB, "w");

#endif
    write2Darray(fichier_Mask, data->coloring, 1, m - 1, 1, n - 1);
    write2Darray(fichier_nx, data->nSx, 1, m-1, 1, n-1);
    write2Darray(fichier_ny, data->nSy, 1, m-1, 1, n-1);

    write2Darray(fichier_U, data->u_n, 0, m - 1, 1, n - 1);
    write2Darray(fichier_V, data->v_n, 1, m - 1, 0, n - 1);
    write2Darray(fichier_Vtx, data->omega, 1, m - 1, 1, n - 1);
    //write2Darray(fichier_P, data->P,1,m-1,1,n-1);
#ifdef TEMP
    write2Darray(fichier_T, data->T_n,1,m-1,1,n-1);
    write2Darray(fichier_YA, data->C_n[0],1,m-1,1,n-1);
    //write2Darray(fichier_YB, data->C_n[1],1,m-1,1,n-1);
#endif

    //CLOSE FILES 
    fclose(fichier_Mask);
    fclose(fichier_nx);
    fclose(fichier_ny);
    fclose(fichier_U);
    fclose(fichier_V);
    fclose(fichier_Vtx);
    //fclose(fichier_P);

#ifdef TEMP
    fclose(fichier_T);
    fclose(fichier_YA);
   // fclose(fichier_YB);
#endif

}

void writeFields_periodic(Data* data, char* folder, int it) {
    int m = data->m;
    int n = data->n;

    FILE* fichier_U = NULL;
    FILE* fichier_V = NULL;
    FILE* fichier_T = NULL;
    FILE* fichier_YA = NULL;
    FILE* fichier_YB = NULL;
    FILE* fichier_Mask = NULL;
    FILE* fichier_Vtx = NULL;
    FILE* fichier_P = NULL;


    char buffer[10];
    sprintf(buffer, "%d", it);

    char fileMask[30];
    strcpy(fileMask, folder);
    strcat(fileMask, "/Mask");
    strcat(fileMask, "-");
    strcat(fileMask, buffer);
    strcat(fileMask, ".txt");
    fichier_Mask = fopen(fileMask, "w");

    char fileU[30];
    strcpy(fileU, folder);
    strcat(fileU, "/U");
    strcat(fileU, "-");
    strcat(fileU, buffer);
    strcat(fileU, ".txt");
    fichier_U = fopen(fileU, "w");

    char fileV[30];
    strcpy(fileV, folder);
    strcat(fileV, "/V");
    strcat(fileV, "-");
    strcat(fileV, buffer);
    strcat(fileV, ".txt");
    fichier_V = fopen(fileV, "w");

    char fileP[30];
    strcpy(fileP, folder);
    strcat(fileP, "/P");
    strcat(fileP, "-");
    strcat(fileP, buffer);
    strcat(fileP, ".txt");
    fichier_P = fopen(fileP, "w");

    char fileVtx[30];
    strcpy(fileVtx, folder);
    strcat(fileVtx, "/Vtx");
    strcat(fileVtx, "-");
    strcat(fileVtx, buffer);
    strcat(fileVtx, ".txt");
    fichier_Vtx = fopen(fileVtx, "w");

#ifdef TEMP
    char fileT[30];
    strcpy(fileT, folder);
    strcat(fileT, "/T");
    strcat(fileT, "-");
    strcat(fileT, buffer);
    strcat(fileT, ".txt");
    fichier_T = fopen(fileT, "w");

    char fileYA[30];
    strcpy(fileYA, folder);
    strcat(fileYA, "/CA");
    strcat(fileYA, "-");
    strcat(fileYA, buffer);
    strcat(fileYA, ".txt");
    fichier_YA = fopen(fileYA, "w");

    char fileYB[30];
    strcpy(fileYB, folder);
    strcat(fileYB, "/CB");
    strcat(fileYB, "-");
    strcat(fileYB, buffer);
    strcat(fileYB, ".txt");
    fichier_YB = fopen(fileYB, "w");
#endif
    write2Darray(fichier_Mask, data->coloring, 0, m, 1, n-1);
    write2Darray(fichier_U, data->u_n, 0, m, 1, n-1);
    write2Darray(fichier_V, data->v_n, 0, m, 0, n-1);
    write2Darray(fichier_P, data->P, 0, m, 1, n-1);
    write2Darray(fichier_Vtx, data->omega, 0, m, 1, n-1);
#ifdef TEMP
    write2Darray(fichier_T, data->T_n,0,m,1,n-1);
    write2Darray(fichier_YA, data->C_n[0],0,m,1,n-1);
    write2Darray(fichier_YB, data->C_n[1],0,m,1,n-1);
#endif

    //CLOSE FILES
    fclose(fichier_Mask);
    fclose(fichier_U);
    fclose(fichier_V);
    fclose(fichier_Vtx);
    fclose(fichier_P);
#ifdef TEMP
    fclose(fichier_T);
    fclose(fichier_YA);
    fclose(fichier_YB);
#endif

}

void write2Darray(FILE* file, double **data, int iStart, int iEnd, int jStart, int jEnd) {
    if (file != NULL){
        for (int i = iStart; i<iEnd; i++){
            int j;
            for (j = jStart; j<jEnd-1; j++) {
                fprintf(file, "%3.13e \t", data[i][j]);
            }
            fprintf(file, "%3.13e", data[i][j]);
            fprintf(file, "\n");
            fflush(file);
        }
    }
    else{
        printf(" An error occured in write2Darray : invalid file pointer.\n");
    }
}

void writeStatistics(Data* data, FILE* file)
{
    fprintf(file, "%3.6e \t  %3.6e \t  %3.6e \n", data->CFL_max, data->Reh_max, data->Reh_omega_max);
    fflush(file);
}
void writeParticle(Data* data, FILE** file_array, int k)
{
    fprintf(file_array[k], "%3.6e \t  %3.6e \t %3.6e \t %3.6e \t %3.6e \t %3.6e \t %3.6e \n",
            data->xg[k][2], data->yg[k][2], data->theta[k][2], data->Up[k][2], data->Vp[k][2], data->Omega_p[k][2], data->Tp[k]);
    fflush(file_array[k]);
}

void writeForces(Data* data, FILE** file_array, int k)
{
    fprintf(file_array[k], "%3.6e \t  %3.6e \t  %3.6e \n", data->Fx[k], data->Fy[k], data->Tz[k]);
    fflush(file_array[k]);
}

void writeFluxes(Data* data, FILE** file_array, int k)
{
    fprintf(file_array[k], "%3.6e \t  %3.6e \n", data->Q[k], data->Qm[k][0]);
    fflush(file_array[k]);
}
