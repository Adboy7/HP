#ifndef MY_HEADER_IMPLICIT  /* Evite la définition multiple */
#define MY_HEADER_IMPLICIT

#include <stdio.h>
#include <stdlib.h>
#include "project2_DONVAL_CHIODO_map.h"
#include "project2_DONVAL_CHIODO_main.h"

void implicit_calcul_eta(double n, double i, double j, struct Param* param, struct Map* map);

#endif