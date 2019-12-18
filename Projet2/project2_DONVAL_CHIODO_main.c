#include "project2_DONVAL_CHIODO_main.h"
#include "project2_DONVAL_CHIODO_map.h"
#include "project2_DONVAL_CHIODO_explicit.h"
#include "project2_DONVAL_CHIODO_implicit.h"
#include <mpi.h>
// #include <omp.h>

struct Param extract_parametres(char* file) {
    FILE *fp;
    struct Param param;

    fp = fopen(file, "r");

    if(fp == NULL) {
        printf("Erreur lors de l'ouverture du fichier des mapètres\n");
        exit(-2);
    }

    fscanf(fp, "%lf", &param.g);
    fscanf(fp, "%lf", &param.gamma);
    fscanf(fp, "%lf", &param.delta_x);
    fscanf(fp, "%lf", &param.delta_y);
    fscanf(fp, "%lf", &param.delta_t);
    fscanf(fp, "%lf", &param.t_max);
    fscanf(fp, "%lf", &param.A);
    fscanf(fp, "%lf", &param.f);
    fscanf(fp, "%u", &param.S);
    fscanf(fp, "%u", &param.s);
    fscanf(fp, "%lf", &param.r_threshold);
    fscanf(fp, "%d", &param.v_condition);
    // fermeture du fichier
    fclose(fp);

    return param;
}

int main(int argc, char* argv[]) {
    if(argc < 4) {
        printf("Arguement(s) manquant\n");
        exit(-1);
    }

    MPI_Init(&argc, &argv);

    // Variables
    struct Map map;
    struct Param param;
    int myrank, np;
    double start_j, end_j, j;
    MPI_Status mystatus1;

    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    param = extract_parametres(argv[1]);

    if(myrank == 0) {
        printf("Les parametres pour la simulation sont les suivants :\n");
        printf("g = %lf\n", param.g);
        printf("y = %lf\n", param.gamma);
        printf("delta_x = %lf\n", param.delta_x);
        printf("delta_y = %lf\n", param.delta_y);
        printf("delta_t = %lf\n", param.delta_t);
        printf("t_max = %lf\n", param.t_max);
        printf("A = %lf\n", param.A);
        printf("f = %lf\n", param.f);
        printf("S = %u\n", param.S);
        printf("s = %u\n", param.s);
        printf("r_threshold = %lf\n", param.r_threshold);
    }

    map = extract_bathymetric(argv[2], param);

    double n = 1.0;
    unsigned int i = 0, index_i, index_j;
    start_j = map.size_u_j * myrank / np + 1;
    end_j = map.size_u_j * (myrank + 1) / np;
    while(n <= (param.t_max / param.delta_t)) {
        if(myrank == 0) {
            printf("n = %lf\n", n);
        }
        j = start_j;
        // if(end_j == map.size_v_j) {
        //     end_j--;
        // }
        // printf("end_j = ")
        while(j < end_j) {
            // printf("j = %lf\n", j);
            i = 0.0;
            // #pragma omp parallel for default(none) private() shared(param, map, j, i, n)
            for(i = 0; i < (unsigned int) map.size_v_i; i++) {
                if(atoi(argv[3]) == 0) {
                    explicit_calcul_eta(n, i, j, &param, &map);
                } else {
                    implicit_calcul_eta(n, i, j, &param, &map);
                }
            }
            j += 1.0;
        }
        n += 1.0;
        
        // writeStepFile(map, n, start_j, end_j);

        int flag_start = 1;
        int flag_ack_fin = 2;

        // SYNCHRO
        if(myrank == 0) {
            printf("debut synchro\n");
            start_j = map.size_u_j * myrank / np + 1;
            end_j = map.size_u_j * (myrank + 1) / np;
            writeStepFile(map, n, start_j, end_j, myrank);
            for(int r = 1; r < np; r++) {
                start_j = map.size_v_j * r / np + 1;
                MPI_Send(&flag_start, 1, MPI_INT, r, 1, MPI_COMM_WORLD);
                if(r < np - 1) {
                    MPI_Recv(&j, 1, MPI_DOUBLE, r, 1, MPI_COMM_WORLD, &mystatus1);
                    MPI_Recv(map.v[convert(j)], map.size_v_i, MPI_DOUBLE, r, 1, MPI_COMM_WORLD, &mystatus1);
                }
                MPI_Send(map.v[convert(start_j)], map.size_v_i, MPI_DOUBLE, r, 1, MPI_COMM_WORLD);
                MPI_Recv(&flag_ack_fin, 1, MPI_INT, r, 1, MPI_COMM_WORLD, &mystatus1);
            }
        } else {

            MPI_Recv(&flag_start, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &mystatus1);
            if(myrank < np - 1) {
                MPI_Send(&end_j, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
                MPI_Send(map.v[convert(end_j)], map.size_v_i, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
            }
            MPI_Recv(map.v[convert(start_j)], map.size_v_i, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &mystatus1);
            writeStepFile(map, n, start_j, end_j, myrank);
            MPI_Send(&flag_ack_fin, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
        }
    }

    destroyMap(map);

    MPI_Finalize();
    
    return 0;
}
