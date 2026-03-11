#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include "utils.h"

// Interpolation (Serial Code)
void interpolation(double *mesh_value, Points *points) {
    double inv_dx = NX;
    double inv_dy = NY;

    for (int p = 0; p < NUM_Points; p++) {
        double x = points->x[p];
        double y = points->y[p];

        int i = (int)(x * inv_dx);
        int j = (int)(y * inv_dy);

        if (i >= NX) i = NX - 1;
        if (j >= NY) j = NY - 1;

        double dx_local = (x - i * dx);
        double dy_local = (y - j * dy);

        double one_minus_dx = dx - dx_local;
        double one_minus_dy = dy - dy_local;

        int base = j * GRID_X + i;

        mesh_value[base]                 += one_minus_dx * one_minus_dy;
        mesh_value[base + 1]             += dx_local * one_minus_dy;
        mesh_value[base + GRID_X]        += one_minus_dx * dy_local;
        mesh_value[base + GRID_X + 1]    += dx_local * dy_local;
    }
}

// Stochastic Mover (Serial Code) 
void mover_serial(Points *points, double deltaX, double deltaY) {
    // Using rand_r to generate random values in serial as well
    unsigned int seed = 123456789;  // Fixed seed for consistency
    for (long long p = 0; p < NUM_Points; p++) {
        double new_x, new_y;

        while (1) {
            double rx = ((double)rand_r(&seed) / RAND_MAX) * 2.0 * deltaX - deltaX;
            double ry = ((double)rand_r(&seed) / RAND_MAX) * 2.0 * deltaY - deltaY;

            new_x = points->x[p] + rx;
            new_y = points->y[p] + ry;

            if (new_x >= 0.0 && new_x <= 1.0 && new_y >= 0.0 && new_y <= 1.0)
                break;
        }

        points->x[p] = new_x;
        points->y[p] = new_y;
    }
}

// Stochastic Mover (Parallel Code) 
void mover_parallel(Points *points, double deltaX, double deltaY) {
    #pragma omp parallel
    {
        // Each thread gets its own seed for random number generation
        unsigned int seed = 123456789 + omp_get_thread_num();  // Better seed initialization

        #pragma omp for schedule(static)
        for (long long p = 0; p < NUM_Points; p++) {
            // Generate displacement that stays inside the domain
            double min_dx = -points->x[p];
            double max_dx = 1.0 - points->x[p];
            double rx = ((double)rand_r(&seed) / RAND_MAX) * (max_dx - min_dx);  // Using rand_r for randomness

            double min_dy = -points->y[p];
            double max_dy = 1.0 - points->y[p];
            double ry = ((double)rand_r(&seed) / RAND_MAX) * (max_dy - min_dy);  // Using rand_r for randomness

            points->x[p] += rx;
            points->y[p] += ry;
        }
    }
}

// Write mesh to file
void save_mesh(double *mesh_value) {
    FILE *fd = fopen("Mesh.out", "w");
    if (!fd) {
        printf("Error creating Mesh.out\n");
        exit(1);
    }

    for (int i = 0; i < GRID_Y; i++) {
        for (int j = 0; j < GRID_X; j++) {
            fprintf(fd, "%lf ", mesh_value[i * GRID_X + j]);
        }
        fprintf(fd, "\n");
    }

    fclose(fd);
}
