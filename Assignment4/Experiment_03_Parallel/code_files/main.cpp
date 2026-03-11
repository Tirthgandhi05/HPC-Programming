#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include "init.h"
#include "utils.h"

// Global variables
int GRID_X, GRID_Y, NX, NY;
int NUM_Points;
int Maxiter;
double dx, dy;

int main() {

    NX = 1000;
    NY = 400;
    NUM_Points = 14000000;
    Maxiter = 10;

    GRID_X = NX + 1;
    GRID_Y = NY + 1;

    dx = 1.0 / NX;
    dy = 1.0 / NY;

    omp_set_num_threads(4); 
    
    double *mesh_value = (double*) calloc(GRID_X * GRID_Y, sizeof(double));
    Points points;
    points.x = (double*) malloc(NUM_Points * sizeof(double));
    points.y = (double*) malloc(NUM_Points * sizeof(double));

    initializepoints(&points);

    printf("Iteration\tInterp(s)\tMover(s)\tTotal(s)\n");

    double serial_mover_total = 0.0;

    // Open the CSV file for writing once at the beginning
    FILE *fp = fopen("exp3_results.csv", "w");
    if (!fp) {
        printf("Error creating CSV file\n");
        exit(1);
    }

    // Write the header for the CSV file
    fprintf(fp, "Iteration,InterpTime,MoverTime,TotalTime\n");

    for (int iter = 0; iter < Maxiter; iter++) {

        // Reset mesh
        for (int i = 0; i < GRID_X * GRID_Y; i++)
            mesh_value[i] = 0.0;

        // Interpolation timing
        double start_interp = omp_get_wtime();
        interpolation(mesh_value, &points);
        double end_interp = omp_get_wtime();
        double interp_time = end_interp - start_interp;

        double start_move = omp_get_wtime();

        if (iter == 0) {
            // Run serial mover once to compute speedup
            mover_serial(&points, dx, dy);
        } else {
            mover_parallel(&points, dx, dy);
        }

        double end_move = omp_get_wtime();
        double move_time = end_move - start_move;

        if (iter == 0) serial_mover_total = move_time;

        double total_time = interp_time + move_time;

        // Print iteration timings to the console
        printf("%d\t\t%lf\t%lf\t%lf\n", iter + 1, interp_time, move_time, total_time);

        // Write the actual times to the CSV file
        fprintf(fp, "%d,%lf,%lf,%lf\n", iter + 1, interp_time, move_time, total_time);
    }

    fclose(fp);  // Close the CSV file after all iterations

    double start_parallel = omp_get_wtime();
    mover_parallel(&points, dx, dy);
    double end_parallel = omp_get_wtime();
    double parallel_time = end_parallel - start_parallel;

    double speedup = serial_mover_total / parallel_time;
    printf("\nSerial Mover Time (1st iteration): %lf s\n", serial_mover_total);
    printf("Parallel Mover Time: %lf s\n", parallel_time);
    printf("Speedup: %lf\n", speedup);

    // Free allocated memory
    free(mesh_value);
    free(points.x);
    free(points.y);

    return 0;
}
