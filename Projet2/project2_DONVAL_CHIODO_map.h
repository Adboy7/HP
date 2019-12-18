#ifndef MY_HEADER_MAP  /* Evite la définition multiple */
#define MY_HEADER_MAP

#include <stdio.h>
#include <stdlib.h>
#include "project2_DONVAL_CHIODO_main.h"


struct Map {
    double a, b;
    unsigned int X, Y;
    double** h;
    double** u;
    double** v;
    double** eta;
    unsigned int size_u_i, size_u_j, size_v_i, size_v_j, size_eta_i, size_eta_j;
};

struct Map extract_bathymetric(char* file, struct Param param);

void destroyMap(struct Map param);

// void printMap(struct Map map);

unsigned int convert(double x);

void print(struct Map map, unsigned int size_i, unsigned int size_j);

void writeStepFile(struct Map map, unsigned int step, double start_j, double end_j, int rank);

double computeSingleH(struct Map *map, struct Param *param, double xw, double yw);

#endif
