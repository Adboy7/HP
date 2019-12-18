#ifndef MY_HEADER_EXPLICIT  /* Evite la définition multiple */
#define MY_HEADER_EXPLICIT

#include <stdio.h>
#include <stdlib.h>
#include "project2_DONVAL_CHIODO_map.h"
#include "project2_DONVAL_CHIODO_main.h"

void explicit_calcul_eta(double n, double i, double j, struct Param* param, struct Map* map);

double explicit_calcul_u(double n, double i, double j, struct Param* param, struct Map* map);

double explicit_calcul_v(double n, double i, double j, struct Param* param, struct Map* map);

#endif