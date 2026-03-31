#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "utils.h"

// -------------------------------------------------------------
// ULTRA-FAST INLINE PRNG (Xorshift32)
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
    return ((double)fast_rand(state) / 4294967295.0) * 2.0 - 1.0;
}

// Returns a double strictly in [0.0, 1.0]
static inline double fast_rand_unit(uint32_t *state) {
    return (double)fast_rand(state) / 4294967295.0;
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

        // Boundary safety
        if (i >= NX) i = NX - 1;
        if (j >= NY) j = NY - 1;

        double dx_local = (x - i * dx);
        double dy_local = (y - j * dy);

        double one_minus_dx = dx - dx_local;
        double one_minus_dy = dy - dy_local;

        int base = j * GRID_X + i;

        // Accumulate weights
        mesh_value[base]                 += one_minus_dx * one_minus_dy;
        mesh_value[base + 1]             += dx_local * one_minus_dy;
        mesh_value[base + GRID_X]        += one_minus_dx * dy_local;
        mesh_value[base + GRID_X + 1]    += dx_local * dy_local;
    }
}

// -------------------------------------------------------------
// APPROACH 1: Immediate Replacement (The Fastest Method)
// -------------------------------------------------------------
void mover_immediate(Points *points, double deltaX, double deltaY) {
    // Static ensures the seed persists across iteration calls
    static uint32_t seed = 123456789; 

    for (int p = 0; p < NUM_Points; p++) {
        double new_x = points->x[p] + fast_rand_norm(&seed) * deltaX;
        double new_y = points->y[p] + fast_rand_norm(&seed) * deltaY;

        // If out of bounds, replace immediately at the same memory index
        if (new_x < 0.0 || new_x > 1.0 || new_y < 0.0 || new_y > 1.0) {
            points->x[p] = fast_rand_unit(&seed);
            points->y[p] = fast_rand_unit(&seed);
        } else {
            points->x[p] = new_x;
            points->y[p] = new_y;
        }
    }
}

// -------------------------------------------------------------
// APPROACH 2: Deferred Insertion (Slower due to Memory Compaction)
// -------------------------------------------------------------
void mover_deferred(Points *points, double deltaX, double deltaY) {
    static uint32_t seed = 987654321; 
    
    int valid_count = 0; // Tracks the write index

    // Pass 1: Move particles and compact valid ones to the front of the array
    for (int p = 0; p < NUM_Points; p++) {
        double new_x = points->x[p] + fast_rand_norm(&seed) * deltaX;
        double new_y = points->y[p] + fast_rand_norm(&seed) * deltaY;

        if (new_x >= 0.0 && new_x <= 1.0 && new_y >= 0.0 && new_y <= 1.0) {
            points->x[valid_count] = new_x;
            points->y[valid_count] = new_y;
            valid_count++;
        }
    }

    // Pass 2: The remaining space (voids) are pushed to the end. Insert new particles here.
    for (int p = valid_count; p < NUM_Points; p++) {
        points->x[p] = fast_rand_unit(&seed);
        points->y[p] = fast_rand_unit(&seed);
    }
}