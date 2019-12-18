#include "project2_DONVAL_CHIODO_explicit.h"
#include <math.h>

void explicit_calcul_eta(double n, double i, double j, struct Param* param, struct Map* map) {
    double p1, p2;
    if(n == 0.0) {
        map->eta[convert(j)][convert(i)] = 0.0;
    } else {
        p1 = (map->h[convert(j)][convert(i + 0.5)] * explicit_calcul_u(n + 0.5, i + 0.5, j, param, map) - map->h[convert(j)][convert(i - 0.5)] * explicit_calcul_u(n + 0.5, i - 0.5, j, param, map)) / param->delta_x;
        p2 = (map->h[convert(j + 0.5)][convert(i)] * explicit_calcul_v(n + 0.5, i, j + 0.5, param, map) - map->h[convert(j - 0.5)][convert(i)] * explicit_calcul_v(n + 0.5, i, j - 0.5, param, map)) / param->delta_y;
        map->eta[convert(j)][convert(i)] = param->delta_t * (-1 * p1 -1 * p2) + map->u[convert(j)][convert(i)] + map->eta[convert(j)][convert(i)];
    }
    // printf("n = %f, eta = %f, i = %d, j = %d\n", n, map->eta[convert(j)][convert(i)], convert(i), convert(j));
}

double explicit_calcul_u(double n, double i, double j, struct Param* param, struct Map* map) {

    if(i == -0.5 || i == map->size_u_i) {
        return 0.0;
    }

    if(n == 0.5) {
        map->u[convert(j)][convert(i + 0.5)] = 0.0;
    } else {
        // printf("pass\n");
        map->u[convert(j)][convert(i + 0.5)] = param->delta_t * (-1 * param->g * ((map->eta[convert(j)][convert(i + 1)] - map->eta[convert(j)][convert(i)])/param->delta_x) -1 * param->gamma * map->u[convert(j)][convert(i)]) + map->u[convert(j)][convert(i + 0.5)];
    }

    // printf("n = %f, u = %f\n", n, map->u[convert(j)][convert(i + 0.5)]);

    return map->u[convert(j)][convert(i + 0.5)];
}

double explicit_calcul_v(double n, double i, double j, struct Param* param, struct Map* map) {

    if(j == -0.5) {
        return 0.0;
    } else if(j == (map->b / param->delta_y) + 0.5) {
        double t = (n + 0.5) / param->delta_t;
        if(param->v_condition == 0) {
            map->v[convert(j)][convert(i)] = param->A * sin(2 * 3.14 * param->f * t);
            return map->v[convert(j)][convert(i)];
        } else {
            map->v[convert(j)][convert(i)] = param->A * sin(2 * 3.14 * param->f * t) * exp(-(t/500));
            return map->v[convert(j)][convert(i)];
        }
    } else if(n == 0.5) {
        map->v[convert(j + 0.5)][convert(i)] = 0.0;
    } else {
        map->v[convert(j + 0.5)][convert(i)] = param->delta_t * (-1 * param->g * ((map->eta[convert(j + 1)][convert(i)] - map->eta[convert(j)][convert(i)])/param->delta_x) -1 * param->gamma * map->v[convert(j)][convert(i)]) + map->v[convert(j + 0.5)][convert(i)];
    }

    return map->v[convert(j + 0.5)][convert(i)];
}