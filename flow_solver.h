//
// Created by Baptiste Hardy on 5/07/18.
//

#ifndef INC_2DINCOMP_FLOW_SOLVER_H
#define INC_2DINCOMP_FLOW_SOLVER_H

void get_masks(Data* data);
void get_Cs(Data* data);
void get_Ts(Data* data);
void get_Us_Vs(Data* data);
void get_Ustar_Vstar(Data* data, double ramp);
void get_vorticity(Data* data);
void get_diffusivity(Data* data);
void get_conductivity(Data* data);
void update_flow(Data* data);
void update_scalars(Data* data);
void get_rate(Data* data, double* r, double*** Cs, double** Ts, int i, int j);
void track_interface(Data* data, int* K, double* THETA, int* right, int* left, int* above, int* below, int i, int j);

#endif //INC_2DINCOMP_FLOW_SOLVER_H
