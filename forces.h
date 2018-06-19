//
// Created by Baptiste Hardy on 8/06/18.
//

#ifndef INC_2DINCOMP_FORCES_H
#define INC_2DINCOMP_FORCES_H

int integrate_penalization(Data *data, double* surf, int k);
void compute_forces_fluxes(Data* data, int k);
void get_tau(Data* data);
void get_tau_periodic(Data* data);
void compute_forces_NOCA(Data* data, FILE* file, int I1, int I2, int J1, int J2);

#endif //INC_2DINCOMP_FORCES_H

