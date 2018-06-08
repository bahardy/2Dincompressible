//
// Created by Baptiste Hardy on 8/06/18.
//

#ifndef INC_2DINCOMP_FORCES_H
#define INC_2DINCOMP_FORCES_H

int integrate_penalization(Data *data, double* surf, int k);
void compute_forces_fluxes(Data* data, int k);

#endif //INC_2DINCOMP_FORCES_H

