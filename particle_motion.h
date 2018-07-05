//
// Created by Baptiste Hardy on 29/06/18.
//

#ifndef INC_2DINCOMP_PARTICLE_MOTION_H
#define INC_2DINCOMP_PARTICLE_MOTION_H

#include "main.h"

void update_Xp(Data* data, int k);
void update_Up(Data* data, int k);
void relax_Up(Data* data, double relax, double* Up_old, double* Vp_old, double* Omega_p_old, int k);
void update_Tp(Data* data,int k);
void update_Cp(Data* data, int k);
void compute_Qr(double** Qr, double rate, double dH, int k);

#endif //INC_2DINCOMP_PARTICLE_MOTION_H
