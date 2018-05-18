//
// Created by Baptiste Hardy on 9/03/18.
//
#include "main.h"
#ifndef INCOMP_POISSON_H
#define INCOMP_POISSON_H

void old_poisson_solver(Data* data);
void old_poisson_solver_periodic(Data* data);
PetscErrorCode poisson_solver(Data* data, int myrank, int nbproc);
PetscErrorCode poisson_solver_periodic(Data* data, int myrank, int nbproc);


#endif //INCOMP_POISSON_H
