//
// Created by Baptiste Hardy on 9/03/18.
//
#include "main.h"
#ifndef 2DINCOMP_POISSON
#define 2DINCOMP_POISSON

void old_poisson_solver(Data* data);
PetscErrorCode poisson_solver(Data* data, int myrank, int nbproc);

#endif //2DINCOMP_POISSON
