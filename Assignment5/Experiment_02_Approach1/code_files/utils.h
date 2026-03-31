#ifndef UTILS_H
#define UTILS_H

#include <stdint.h>

typedef struct {
    double *x;
    double *y;
} Points;

// Fast RNG
uint32_t fast_rand(uint32_t *state);
double fast_rand_unit(uint32_t *state);
double fast_rand_norm(uint32_t *state);

// Interpolation (Optimized without Atomics)
void interpolation_parallel(double *mesh, Points *points, int NX, int NY, int GRID_X, int N, double dx, double dy, int threads);

// Serial Baselines (For accurate speedup calculation)
void mover_immediate_serial(Points *points, int N, double dx, double dy);
void mover_deferred_serial(Points *points, int N, double dx, double dy);
void mover_baseline_serial(Points *points, int N, double dx, double dy);

// Parallel Movers
void mover_immediate_parallel(Points *points, int N, double dx, double dy, int threads);
void mover_deferred_parallel(Points *points, int N, double dx, double dy, int threads);
void mover_baseline_parallel(Points *points, int N, double dx, double dy, int threads);

#endif