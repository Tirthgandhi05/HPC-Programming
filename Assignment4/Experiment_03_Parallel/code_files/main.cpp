#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include "init.h"
#include "utils.h"

int GRID_X, GRID_Y, NX, NY;
int NUM_Points;
int Maxiter;
double dx, dy;

int main() {
    NX = 1000;
    NY = 400;
    NUM_Points = 14000000; // 14 Million Particles
    Maxiter = 10;

    GRID_X = NX + 1;
    GRID_Y = NY + 1;
    dx = 1.0 / NX;
    dy = 1.0 / NY;

    // Fixed to 4 threads as per manual
    omp_set_num_threads(4);
   
    double *mesh_value = (double*) calloc(GRID_X * GRID_Y, sizeof(double));
    Points points;
    points.x = (double*) malloc(NUM_Points * sizeof(double));
    points.y = (double*) malloc(NUM_Points * sizeof(double));

    initializepoints(&points);

    // ========================================================
    // Phase 1: Get an accurate Serial Mover baseline
    // ========================================================
    printf("Running Serial Mover baseline...\n");
    double start_serial = omp_get_wtime();
    for (int i = 0; i < Maxiter; i++) {
        mover_serial(&points, dx, dy);
    }
    double end_serial = omp_get_wtime();
    double avg_serial_time = (end_serial - start_serial) / Maxiter;


    // ========================================================
    // Phase 2: Run the 10 Iterations (Interp + Parallel Mover)
    // ========================================================
    FILE *fp = fopen("exp3_results.csv", "w");
    fprintf(fp, "Iteration,InterpTime,MoverTime,TotalTime\n");
   
    printf("\nIteration\tInterp(s)\tMover_Par(s)\tTotal(s)\n");

    double total_parallel_mover_time = 0.0;

    for (int iter = 0; iter < Maxiter; iter++) {
        for (int i = 0; i < GRID_X * GRID_Y; i++) mesh_value[i] = 0.0;

        // 1. Interpolation
        double start_interp = omp_get_wtime();
        interpolation(mesh_value, &points);
        double end_interp = omp_get_wtime();
        double interp_time = end_interp - start_interp;

        // 2. Parallel Mover
        double start_move = omp_get_wtime();
        mover_parallel(&points, dx, dy);
        double end_move = omp_get_wtime();
        double move_time = end_move - start_move;
       
        total_parallel_mover_time += move_time;
        double total_time = interp_time + move_time;

        printf("%d\t\t%.4lf\t\t%.4lf\t\t%.4lf\n", iter + 1, interp_time, move_time, total_time);
        fprintf(fp, "%d,%lf,%lf,%lf\n", iter + 1, interp_time, move_time, total_time);
    }
    fclose(fp);

    double avg_parallel_time = total_parallel_mover_time / Maxiter;
    double speedup = avg_serial_time / avg_parallel_time;

    // Save speedup results for python script
    FILE *fp2 = fopen("speedup_data.txt", "w");
    fprintf(fp2, "%lf\n%lf\n%lf\n", avg_serial_time, avg_parallel_time, speedup);
    fclose(fp2);

    printf("\n--- Final Results ---\n");
    printf("Avg Serial Mover Time   : %.4lf s\n", avg_serial_time);
    printf("Avg Parallel Mover Time : %.4lf s\n", avg_parallel_time);
    printf("Achieved Speedup (S)    : %.4lf x\n", speedup);

    free(mesh_value);
    free(points.x);
    free(points.y);
    return 0;
}






