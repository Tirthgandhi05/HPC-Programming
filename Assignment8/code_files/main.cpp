#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "init.h"
#include "utils.h"

/* Global variable definitions (declared extern in init.h) */
int    GRID_X, GRID_Y, NX, NY;
int    NUM_Points, Maxiter;
double dx, dy;

/* Wall-clock timer — matches MPI_Wtime() semantics */
static inline double now_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

int main(int argc, char **argv)
{
    if (argc < 2 || argc > 3) {
        printf("Usage: %s <input_file>\n", argv[0]);
        return 1;
    }

    FILE *file = fopen(argv[1], "rb");
    if (!file) {
        printf("Error: cannot open '%s'\n", argv[1]);
        return 1;
    }

    /* Read header */
    fread(&NX,         sizeof(int), 1, file);
    fread(&NY,         sizeof(int), 1, file);
    fread(&NUM_Points, sizeof(int), 1, file);
    fread(&Maxiter,    sizeof(int), 1, file);

    GRID_X = NX + 1;
    GRID_Y = NY + 1;
    dx     = 1.0 / NX;
    dy     = 1.0 / NY;

    printf("Serial PIC | NX=%d NY=%d Points=%d Maxiter=%d\n",
           NX, NY, NUM_Points, Maxiter);

    double *mesh_value = (double *)calloc(GRID_X * GRID_Y, sizeof(double));
    Points *points     = (Points *)malloc(NUM_Points * sizeof(Points));
    if (!mesh_value || !points) { printf("Alloc failed\n"); return 1; }

    /* Read initial positions once, before the loop */
    read_points(file, points, NUM_Points);
    fclose(file);

    double total_interp = 0.0, total_mover = 0.0;

    int iter;
    for (iter = 0; iter < Maxiter; iter++) {

        /* 1. Forward interpolation */
        memset(mesh_value, 0, GRID_X * GRID_Y * sizeof(double));
        double t0 = now_sec();
        interpolation(mesh_value, points, NUM_Points);
        double t1 = now_sec();
        total_interp += t1 - t0;

        /* 2. Normalise */
        double mesh_min, mesh_max;
        normalise_mesh(mesh_value, &mesh_min, &mesh_max);

        /* 3. Mover */
        double t2 = now_sec();
        mover(mesh_value, points, NUM_Points);
        double t3 = now_sec();
        total_mover += t3 - t2;

        /* 4. Denormalise */
        denormalise_mesh(mesh_value, mesh_min, mesh_max);
    }

    save_mesh(mesh_value);

    printf("Interpolation time = %.6lf s\n", total_interp);
    printf("Mover time         = %.6lf s\n", total_mover);
    printf("Total pipeline time      : %.6lf s\n", total_interp + total_mover);

    free(mesh_value);
    free(points);
    return 0;
}