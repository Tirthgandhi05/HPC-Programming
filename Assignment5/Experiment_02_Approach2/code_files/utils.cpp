#include <stdint.h>
#include <stdlib.h>
#include <omp.h>
#include "utils.h"

// ---------------- RNG ----------------
uint32_t fast_rand(uint32_t *state)
{
    uint32_t x = *state;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    *state = x;
    return x;
}

double fast_rand_unit(uint32_t *state)
{
    return (double)fast_rand(state) / 4294967295.0;
}

double fast_rand_norm(uint32_t *state)
{
    return ((double)fast_rand(state) / 4294967295.0) * 2.0 - 1.0;
}

// ---------------- INTERPOLATION ----------------
void interpolation_parallel(double *mesh, Points *points, int NX, int NY, int GRID_X, int N, double dx, double dy, int threads)
{

    int grid_size = GRID_X * (NY + 1);
    double **local = (double **)malloc(threads * sizeof(double *));

#pragma omp parallel num_threads(threads)
    {
        int tid = omp_get_thread_num();
        local[tid] = (double *)calloc(grid_size, sizeof(double));

#pragma omp for schedule(static)
        for (int p = 0; p < N; p++)
        {
            double x = points->x[p];
            double y = points->y[p];

            int i = (int)(x * NX);
            int j = (int)(y * NY);

            // ✅ FIXED boundary
            if (i >= NX - 1)
                i = NX - 2;
            if (j >= NY - 1)
                j = NY - 2;

            double dx_local = x - i * dx;
            double dy_local = y - j * dy;

            double w1 = (dx - dx_local) * (dy - dy_local);
            double w2 = dx_local * (dy - dy_local);
            double w3 = (dx - dx_local) * dy_local;
            double w4 = dx_local * dy_local;

            int base = j * GRID_X + i;

            local[tid][base] += w1;
            local[tid][base + 1] += w2;
            local[tid][base + GRID_X] += w3;
            local[tid][base + GRID_X + 1] += w4;
        }
    }

    // reduction
    for (int t = 0; t < threads; t++)
    {
        for (int i = 0; i < grid_size; i++)
        {
            mesh[i] += local[t][i];
        }
        free(local[t]);
    }
    free(local);
}

// ---------------- SERIAL BASELINES ----------------
void mover_immediate_serial(Points *points, int N, double dx, double dy)
{
    uint32_t seed = 1234;
    for (int p = 0; p < N; p++)
    {
        double nx = points->x[p] + fast_rand_norm(&seed) * dx;
        double ny = points->y[p] + fast_rand_norm(&seed) * dy;

        if (nx < 0 || nx > 1 || ny < 0 || ny > 1)
        {
            points->x[p] = fast_rand_unit(&seed);
            points->y[p] = fast_rand_unit(&seed);
        }
        else
        {
            points->x[p] = nx;
            points->y[p] = ny;
        }
    }
}

void mover_deferred_serial(Points *points, int N, double dx, double dy)
{
    uint32_t seed = 999;
    int valid = 0;

    for (int p = 0; p < N; p++)
    {
        double nx = points->x[p] + fast_rand_norm(&seed) * dx;
        double ny = points->y[p] + fast_rand_norm(&seed) * dy;

        if (nx >= 0 && nx <= 1 && ny >= 0 && ny <= 1)
        {
            points->x[valid] = nx;
            points->y[valid] = ny;
            valid++;
        }
    }

    for (int p = valid; p < N; p++)
    {
        points->x[p] = fast_rand_unit(&seed);
        points->y[p] = fast_rand_unit(&seed);
    }
}

// ---------------- PARALLEL ----------------
void mover_immediate_parallel(Points *points, int N, double dx, double dy, int threads)
{
#pragma omp parallel num_threads(threads)
    {
        uint32_t seed = 1234 + omp_get_thread_num();

#pragma omp for schedule(static)
        for (int p = 0; p < N; p++)
        {
            double nx = points->x[p] + fast_rand_norm(&seed) * dx;
            double ny = points->y[p] + fast_rand_norm(&seed) * dy;

            if (nx < 0 || nx > 1 || ny < 0 || ny > 1)
            {
                points->x[p] = fast_rand_unit(&seed);
                points->y[p] = fast_rand_unit(&seed);
            }
            else
            {
                points->x[p] = nx;
                points->y[p] = ny;
            }
        }
    }
}

void mover_deferred_parallel(Points *points, int N, double dx, double dy, int threads)
{

    double *temp_x = (double *)malloc(N * sizeof(double));
    double *temp_y = (double *)malloc(N * sizeof(double));
    int *offset = (int *)calloc(threads + 1, sizeof(int));

// Pass 1: count
#pragma omp parallel num_threads(threads)
    {
        int tid = omp_get_thread_num();
        uint32_t seed = 999 + tid;
        int count = 0;

#pragma omp for schedule(static)
        for (int p = 0; p < N; p++)
        {
            double nx = points->x[p] + fast_rand_norm(&seed) * dx;
            double ny = points->y[p] + fast_rand_norm(&seed) * dy;

            if (nx >= 0 && nx <= 1 && ny >= 0 && ny <= 1)
                count++;
        }

        offset[tid + 1] = count;
    }

    for (int i = 0; i < threads; i++)
        offset[i + 1] += offset[i];

// Pass 2: write
#pragma omp parallel num_threads(threads)
    {
        int tid = omp_get_thread_num();
        uint32_t seed = 999 + tid;
        int idx = offset[tid];

#pragma omp for schedule(static)
        for (int p = 0; p < N; p++)
        {
            double nx = points->x[p] + fast_rand_norm(&seed) * dx;
            double ny = points->y[p] + fast_rand_norm(&seed) * dy;

            if (nx >= 0 && nx <= 1 && ny >= 0 && ny <= 1)
            {
                temp_x[idx] = nx;
                temp_y[idx] = ny;
                idx++;
            }
        }
    }

    int valid = offset[threads];

// copy valid
#pragma omp parallel for num_threads(threads)
    for (int i = 0; i < valid; i++)
    {
        points->x[i] = temp_x[i];
        points->y[i] = temp_y[i];
    }

// ✅ FIXED: PARALLEL insertion
#pragma omp parallel for num_threads(threads)
    for (int i = valid; i < N; i++)
    {
        uint32_t seed = 5555 + i;
        points->x[i] = fast_rand_unit(&seed);
        points->y[i] = fast_rand_unit(&seed);
    }

    free(temp_x);
    free(temp_y);
    free(offset);
}