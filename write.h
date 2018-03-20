//
// Created by Baptiste Hardy on 8/03/18.
//

#ifndef WRITE_H
#define WRITE_H

void writeFields(Data* data, int it);
void writeStatistics(Data* data, FILE* file);
void write2Darray(FILE* file, double **data, int iStart, int iEnd, int jStart, int jEnd);
void writeData(FILE* fichier_data, Data data);


#endif //INC_WRITE_H
