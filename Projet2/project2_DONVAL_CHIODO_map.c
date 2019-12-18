#include "project2_DONVAL_CHIODO_map.h"
struct Map extract_bathymetric(char* file, struct Param param) {
    FILE *fp;
    struct Map map;
    double tmp;
    unsigned int i, j;


    fp = fopen(file, "rb");

    if(fp == NULL) {
        printf("Erreur lors de l'ouverture du fichier des mapètres\n");
        exit(-2);
    }

    fread(&map, sizeof(struct Map), 1, fp);
    printf("a = %lf\tb = %lf\tX = %u\tY = %u\n", map.a, map.b, map.X, map.Y);

    // alloc h table
    map.h = (double **)malloc(sizeof(double*) * map.Y);
    if(map.h == NULL) {
        printf("Erreur matrice = NULL\n");
        exit(-4);
    }
    for(i = 0; i < map.Y; i++) {
        // printf("i = %u\n", i);
        map.h[i] = (double *)malloc(sizeof(double) * map.X);
    }

    for(i = 0; i < map.Y; i++) {
        for(unsigned int j = 0; j < map.X; j++) {
            fread(&tmp, sizeof(double), 1, fp);
            // printf("x = %u, y = %u, value = %lf\n", j, i, tmp);

            map.h[i][j] = tmp;
        }
    }

    // fermeture du fichier
    fclose(fp);

    // alloc u table
    map.size_u_j = ((map.b / param.delta_y) -1) * 2;
    map.size_u_i = ((map.a / param.delta_x) -0.5) * 2;
    // printf("size_u_i = %d, size_u_j = %d\n", map.size_u_i, map.size_u_j);
    map.u = (double**)malloc(sizeof(double*) * map.size_u_j);
    for(j = 0; j < map.size_u_j; j++) {
        map.u[j] = (double *)malloc(sizeof(double) * map.size_u_i);
    }

    // alloc v table
    map.size_v_j = ((map.b / param.delta_y) + 1.5) * 2;
    map.size_v_i = ((map.a / param.delta_x) -1) * 2;
    // printf("size_v_i = %d, size_v_j = %d\n", map.size_v_i, map.size_v_j);
    map.v = (double**)malloc(sizeof(double*) * map.size_v_j);
    for(j = 0; j < map.size_v_j; j++) {
        map.v[j] = (double *)malloc(sizeof(double) * map.size_v_i);
    }

    map.eta = (double**)malloc(sizeof(double*) * map.size_v_j);
    for(j = 0; j < map.size_v_j; j++) {
        map.eta[j] = (double *)malloc(sizeof(double) * map.size_u_i);
    }

    double** tmp_h = (double**) malloc(sizeof(double*) * map.size_u_j);

    printf("compute new h\n");
    for(j = 0; j < map.size_u_j; j++) {
        tmp_h[j] = (double*) malloc(sizeof(double) * map.size_v_i);
        for(i = 0; i < map.size_v_i; i++) {
            tmp_h[j][i] = computeSingleH(&map, &param, i, j);
        }
    }

    // printf("compute new h after\n");

    // free old 
    for(j = 0; j < map.Y; j++) {
        free(map.h[j]);
    }
    free(map.h);

    map.h = tmp_h;

    return map;
}

void destroyMap(struct Map map) {
    unsigned int i;
    for(i = 0; i < map.Y; ++i) {
        free(map.h[i]);
    }

    for(i = 0; i < map.size_u_j; i++) {
        free(map.u[i]);
    }

    for(i = 0; i < map.size_v_j; i++) {
        free(map.v[i]);
        free(map.eta[i]);
    }

    free(map.h);
    free(map.u);
    free(map.v);
    free(map.eta);
}

// void printMap(struct Map map) {
//     for(unsigned int i = 0; i < map.Y; i++) {
//         for(unsigned int j = 0; j < map.X; j++) {
//             printf("%lf ", map.h[i][j]);
//         }
//         printf("\n");
//     }
// }

unsigned int convert(double x) {
    if(x == -0.5) {
        return 0;
    } else {
        return (unsigned int) x + 0.5;
    }
}

void print(struct Map map, unsigned int size_i, unsigned int size_j) {
    for(unsigned j = 0; j < size_j; j++) {
        for(unsigned i = 0; i < size_i; i++) {
            printf("%f ", map.eta[j][i]);
        }
        printf("\n");
    }
}

void writeStepFile(struct Map map, unsigned int step, double start_j, double end_j, int rank) {

    char file[100];
    FILE *fp;
    unsigned int i, j;


    unsigned int sizeui = map.size_u_i - 1;
    unsigned int sizeuj = map.size_u_j - 1;
    unsigned int sizevi = map.size_v_i - 1;
    unsigned int sizevj = map.size_v_j - 1;



    // ETA
    sprintf(file, "./results/eta_%u.dat", step);
    fp = fopen(file, "wb");
    
    if(fp == NULL) {
        printf("erreur open %s\n", file);
        return;
    }
    if(rank == 0){
	printf("%u\n", (unsigned int)map.size_v_i);
        //fprintf(fp, "%u ", (unsigned int)map.size_v_i);
        //fprintf(fp, "%u ", (unsigned int)map.size_u_j);
        fwrite(&sizevi, sizeof(unsigned int), 1, fp);
	    fwrite(&sizevj, sizeof(unsigned int), 1, fp);
	    printf("%u,%u pass\n",sizevi,sizevj);
    }
    for(j = convert(start_j); j < convert(end_j); j++) {
        // for(i = 0; i < map.size_v_i; i++) {
            //fprintf(fp,"%lf ",map.eta[j][i]);
            fwrite(map.eta[j], sizeof(double), map.size_v_i, fp);
        // }
    }

    fclose(fp);

    // U
    sprintf(file, "./results/u_%u.dat", step);
    fp = fopen(file, "wb");
    if(fp == NULL) {
        printf("erreur open %s\n", file);
        return;
    }
    if(rank == 0){
        // fprintf(fp, "%u ", (unsigned int)map.size_u_i);
        // fprintf(fp, "%u ", (unsigned int)map.size_u_j);
        fwrite(&sizeui, sizeof(unsigned int), 1, fp);
        fwrite(&sizeuj, sizeof(unsigned int), 1, fp);
    }
    for(j = convert(start_j); j < convert(end_j); j++) {
        // for(i = 0; i < map.size_u_i; i++) {
            //fprintf(fp,"%lf ",map.u[j][i]);
            fwrite(map.u[j], sizeof(double), map.size_u_i, fp);
        // }
    }

    fclose(fp);
    
    // V
    sprintf(file, "./results/v_%u.dat", step);
    fp = fopen(file, "wb");

    if(fp == NULL) {
        printf("erreur open %s\n", file);
        return;
    }
    if(rank == 0){
        // :fprintf(fp, "%u ", (unsigned int)map.size_v_i);
        // fprintf(fp, "%u ", (unsigned int)map.size_v_j);
        fwrite(&sizevi, sizeof(unsigned int), 1, fp);
        fwrite(&sizevj, sizeof(unsigned int), 1, fp);
    }
    for(j = convert(start_j); j < convert(end_j); j++) {
        // for(i = 0; i < map.size_v_i; i++) {
            //fprintf(fp,"%lf ",map.v[j][i]);
            fwrite(map.v[j], sizeof(double), map.size_v_i, fp);
        // }
    }

    fclose(fp);
}

double computeSingleH(struct Map *map, struct Param *param, double xw, double yw)
{
    double stepX = map->a / (map->X - 1);
    double stepY = map->b / (map->Y - 1);
    double xk = (-param->delta_x) / 2;
    double yl = (-param->delta_y) / 2;
    int k = 0;
    double newH;
    while (xk < map->a - (param->delta_x) / 2)
    {
        xk = xk + k * stepX;
        if (xk <= xw && xk + stepX >= xw)
        {
            break;
        }
        k++;
    }
    int l = 0;
    while (yl < map->b - (param->delta_x) / 2)
    {
        yl = yl + l * stepY;
        if (yl <= yw && yl + stepY >= yw)
        {
            break;
        }
        l++;
    }

    newH = ((xk + stepX - xw) * (yl + stepY - yw) * map->h[k][l] + (xk + stepX - xw) * (yw - yl) * map->h[k][l + 1] + (xw - xk) * (yl + stepY - yw) * map->h[k + 1][l] 
        + (xw - xk) * (yw - yl) * map->h[k + 1][l + 1]) / (stepX * stepY);
    
    return newH;
}
