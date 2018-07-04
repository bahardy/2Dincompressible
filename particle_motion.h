//
// Created by Baptiste Hardy on 29/06/18.
//

#ifndef INC_2DINCOMP_PARTICLE_MOTION_H
#define INC_2DINCOMP_PARTICLE_MOTION_H

#include "main.h"


void update_Xp(Data* data, double* Xp_k, double* Yp_k, double* theta_k,
               double* Up_k, double* Vp_k, double* Omega_p_k,  int k);
void update_Up(Data* data, double* Up_k, double* Vp_k, double* Omega_p_k, int k);
#endif //INC_2DINCOMP_PARTICLE_MOTION_H
