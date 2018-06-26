//
// Created by Baptiste Hardy on 8/03/18.
//

#ifndef WRITE_H
#define WRITE_H

#include "main.h"


void writeFields(Data* data, int it);
void writeFields_periodic(Data* data, int i);
void writeMask(Data* data);
void write2Darray(FILE* file, double **data, int iStart, int iEnd, int jStart, int jEnd);
void writeData(FILE* fichier_data, Data data);
void writeStatistics(Data* data, FILE* file);
void writeParticle(Data* data, FILE** file_array, int k);
void writeForces(Data* data, FILE** file_array, int k);
void writeFluxes(Data* data, FILE** file_array, int k);


#endif //INC_WRITE_H
