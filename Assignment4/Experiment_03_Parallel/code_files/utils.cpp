#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <omp.h>
#include "utils.h"
#include "init.h"

// -------------------------------------------------------------
// ULTRA-FAST INLINE PRNG (Xorshift32)
// This avoids the massive CPU overhead of rand() and rand_r()
// -------------------------------------------------------------
static inline uint32_t fast_rand(uint32_t *state) {
    uint32_t x = *state;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    *state = x;
    return x;
}

// Returns a double strictly between -1.0 and 1.0
static inline double fast_rand_norm(uint32_t *state) {
    // 4294967295.0 is the max value of uint32_t
    return ((double)fast_rand(state) / 4294967295.0) * 2.0 - 1.0;
}

// -------------------------------------------------------------
// INTERPOLATION (Serial)
// -------------------------------------------------------------
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

// -------------------------------------------------------------
// STOCHASTIC MOVER (Serial Baseline)
// -------------------------------------------------------------
void mover_serial(Points *points, double deltaX, double deltaY) {
    uint32_t seed = 123456789; // Fixed seed
   
    for (int p = 0; p < NUM_Points; p++) {
        double orig_x = points->x[p];
        double orig_y = points->y[p];
        double new_x, new_y;

        // Rejection Sampling: Keep rolling until it stays inside the domain
        while (1) {
            new_x = orig_x + fast_rand_norm(&seed) * deltaX;
            new_y = orig_y + fast_rand_norm(&seed) * deltaY;
           
            if (new_x >= 0.0 && new_x <= 1.0 && new_y >= 0.0 && new_y <= 1.0) {
                break;
            }
        }

        points->x[p] = new_x;
        points->y[p] = new_y;
    }
}

// -------------------------------------------------------------
// STOCHASTIC MOVER (OpenMP Parallel)
// -------------------------------------------------------------
void mover_parallel(Points *points, double deltaX, double deltaY) {
   
    #pragma omp parallel
    {
        // Give every thread its own unique seed to prevent false sharing/locking
        uint32_t seed = 123456789 + omp_get_thread_num();

        // Static schedule is best here since workload per particle is relatively uniform
        #pragma omp for schedule(static)
        for (int p = 0; p < NUM_Points; p++) {
            double orig_x = points->x[p];
            double orig_y = points->y[p];
            double new_x, new_y;

            while (1) {
                new_x = orig_x + fast_rand_norm(&seed) * deltaX;
                new_y = orig_y + fast_rand_norm(&seed) * deltaY;
               
                if (new_x >= 0.0 && new_x <= 1.0 && new_y >= 0.0 && new_y <= 1.0) {
                    break;
                }
            }

            points->x[p] = new_x;
            points->y[p] = new_y;
        }
    }
}

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

