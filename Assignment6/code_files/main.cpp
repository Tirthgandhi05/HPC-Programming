#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#include "init.h"
#include "utils.h"

int GRID_X, GRID_Y, NX, NY;
int NUM_Points, Maxiter;
double dx, dy;

int main(int argc, char **argv) {
    if (argc != 3) {
        printf("Usage: %s <input_file> <config_name>\n", argv[0]);
        return 1;
    }

    const char* config_name = argv[2];
    FILE *file = fopen(argv[1], "rb");
    if (!file) exit(1);

    fread(&NX, sizeof(int), 1, file);
    fread(&NY, sizeof(int), 1, file);
    fread(&NUM_Points, sizeof(int), 1, file);
    fread(&Maxiter, sizeof(int), 1, file);

    GRID_X = NX + 1;
    GRID_Y = NY + 1;
    dx = 1.0 / NX;
    dy = 1.0 / NY;

    long data_start_offset = ftell(file);

    double *mesh_value = (double *) calloc(GRID_X * GRID_Y, sizeof(double));
    Points *points = (Points *) calloc(NUM_Points, sizeof(Points));
    double *all_local_meshes = (double *)calloc(16 * GRID_X * GRID_Y, sizeof(double));

    FILE *csv = fopen("results.csv", "a");
    if (!csv) exit(1);

    fseek(csv, 0, SEEK_END);
    if (ftell(csv) == 0)
        fprintf(csv, "Config,Threads,Time,Speedup,Efficiency\n");

    int thread_counts[] = {1, 2, 4, 8, 16};
    double serial_time = 0.0;

    for (int t = 0; t < 5; t++) {
        int num_threads = thread_counts[t];
        omp_set_num_threads(num_threads);

        fseek(file, data_start_offset, SEEK_SET);

        double total_time = 0.0;

        for (int iter = 0; iter < Maxiter; iter++) {
            read_points(file, points);

            // ✅ Always reset mesh before each iteration
            memset(mesh_value, 0, GRID_X * GRID_Y * sizeof(double));

            double start = omp_get_wtime();

            if (num_threads == 1) {
                interpolation(mesh_value, points);
            } else {
                parallel_interpolation(mesh_value, points, all_local_meshes, num_threads);
            }

            total_time += (omp_get_wtime() - start);
        }

        double speedup = 1.0, efficiency = 1.0;

        if (num_threads == 1) {
            serial_time = total_time;
            save_mesh(mesh_value); 
        } else {
            speedup = serial_time / total_time;
            efficiency = speedup / num_threads;
        }

        printf("Threads: %2d | Time: %lf s | Speedup: %.2fx | Eff: %.2f\n",
               num_threads, total_time, speedup, efficiency);

        fprintf(csv, "%s,%d,%lf,%lf,%lf\n",
                config_name, num_threads, total_time, speedup, efficiency);
    }

    free(all_local_meshes);
    free(mesh_value);
    free(points);
    fclose(file);
    fclose(csv);

    return 0;
}
