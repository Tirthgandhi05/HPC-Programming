#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "init.h"
#include "utils.h"

int GRID_X, GRID_Y, NX, NY;
int NUM_Points;
int Maxiter;
double dx, dy;

// Particle scaling configurations
int particle_range[5] = {100, 10000, 1000000, 100000000, 1000000000};

void run_experiment(int nx, int ny, int approach)
{
    NX = nx;
    NY = ny;
    GRID_X = NX + 1;
    GRID_Y = NY + 1;
    dx = 1.0 / NX;
    dy = 1.0 / NY;

    char filename[100];
    if (approach == 1)
    {
        sprintf(filename, "exp01_immediate_%dx%d.csv", NX, NY);
        printf("\n--- Running IMMEDIATE Replacement (%dx%d) ---\n", NX, NY);
    }
    else
    {
        sprintf(filename, "exp01_deferred_%dx%d.csv", NX, NY);
        printf("\n--- Running DEFERRED Insertion (%dx%d) ---\n", NX, NY);
    }

    FILE *fp = fopen(filename, "w");
    if (!fp)
    {
        printf("Error: Could not create output file.\n");
        exit(1);
    }
    fprintf(fp, "Particles,InterpTime,MoverTime,TotalTime\n");

    for (int p_idx = 0; p_idx < 5; p_idx++)
    {
        NUM_Points = particle_range[p_idx];

        // Allocate Memory
        double *mesh_value = (double *)calloc(GRID_X * GRID_Y, sizeof(double));
        Points points;
        points.x = (double *)malloc(NUM_Points * sizeof(double));
        points.y = (double *)malloc(NUM_Points * sizeof(double));

        if (!points.x || !points.y || !mesh_value)
        {
            printf("Memory allocation failed at %d particles. Consider reducing limit.\n", NUM_Points);
            break;
        }

        // Initialize points exactly ONCE outside the loop
        initializepoints(&points);

        double interp_total = 0.0;
        double mover_total = 0.0;

        for (int iter = 0; iter < Maxiter; iter++)
        {
            // Zero the mesh safely
            for (int i = 0; i < GRID_X * GRID_Y; i++)
                mesh_value[i] = 0.0;

            double t1 = omp_get_wtime();
            interpolation(mesh_value, &points);
            double t2 = omp_get_wtime();

            if (approach == 1)
            {
                mover_immediate(&points, dx, dy);
            }
            else
            {
                mover_deferred(&points, dx, dy);
            }
            double t3 = omp_get_wtime();

            interp_total += (t2 - t1);
            mover_total += (t3 - t2);
        }

        double total = interp_total + mover_total;

        printf("Particles: %10d | Interp: %.4lf s | Mover: %.4lf s | Total: %.4lf s\n",
               NUM_Points, interp_total, mover_total, total);

        fprintf(fp, "%d,%lf,%lf,%lf\n", NUM_Points, interp_total, mover_total, total);

        free(mesh_value);
        free(points.x);
        free(points.y);
    }
    fclose(fp);
}

int main()
{
    Maxiter = 10;
    int grids[3][2] = {{250, 100}, {500, 200}, {1000, 400}};

    for (int g = 0; g < 3; g++)
    {
        int nx = grids[g][0];
        int ny = grids[g][1];

        // Run Approach 1: Immediate
        run_experiment(nx, ny, 1);
    }

    return 0;
}