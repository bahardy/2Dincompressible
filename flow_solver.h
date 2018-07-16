//
// Created by Baptiste Hardy on 5/07/18.
//

#ifndef INC_2DINCOMP_FLOW_SOLVER_H
#define INC_2DINCOMP_FLOW_SOLVER_H

void get_ghosts(Data* data, double T0, double* C0);
void get_masks(Data* data);
void get_Cs(Data* data);
void get_Ts(Data* data);
void get_Us_Vs(Data* data);
void get_Ustar_Vstar(Data* data, double ramp);
void get_vorticity(Data* data);
void update_flow(Data* data);
void update_scalars(Data* data);

#endif //INC_2DINCOMP_FLOW_SOLVER_H
