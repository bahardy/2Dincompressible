//
// Created by Baptiste Hardy on 4/07/18.
//

#ifndef INC_2DINCOMP_SET_UP_H
#define INC_2DINCOMP_SET_UP_H

#include "main.h"

void set_up(Data *data, int argc, char **argv, int rank);
void set_up_periodic(Data *data, int argc, char **argv, int rank);
void initialize_fields(Data* data);
void initialize_fields_periodic(Data* data);


#endif //INC_2DINCOMP_SET_UP_H
