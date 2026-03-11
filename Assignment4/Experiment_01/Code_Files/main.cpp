#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include "init.h"
#include "utils.h"

int GRID_X, GRID_Y, NX, NY;
int NUM_Points, Maxiter;
double dx, dy;

int main() {

    Maxiter = 10;

    // Particle counts: 10^2 to 10^9
    long long particle_counts[5] = {
	    1LL * 100,
	    1LL * 10000,
	    1LL * 1000000,
	    1LL * 100000000
	};

    // Grid configurations
    int NX_vals[3] = {250, 500, 1000};
    int NY_vals[3] = {100, 200, 400};

    FILE *fp = fopen("exp1_results.csv", "w");
    fprintf(fp, "Config,NX,NY,Particles,TotalInterpolationTime\n");

    for (int config = 0; config < 3; config++) {

        NX = NX_vals[config];
        NY = NY_vals[config];

        GRID_X = NX + 1;
        GRID_Y = NY + 1;

        dx = 1.0 / NX;
        dy = 1.0 / NY;

        printf("\nRunning Configuration %d (NX=%d, NY=%d)\n",
               config+1, NX, NY);

        for (int pcase = 0; pcase < 4; pcase++) {

            NUM_Points = particle_counts[pcase];

            printf("Particles = %lld\n", (long long)NUM_Points);

            // Allocate memory OUTSIDE timing
            double *mesh_value =
                (double*) calloc(GRID_X * GRID_Y, sizeof(double));

            Points points;
            points.x =
                (double*) malloc(NUM_Points * sizeof(double));
            points.y =
                (double*) malloc(NUM_Points * sizeof(double));

            double total_interp_time = 0.0;

            for (int iter = 0; iter < Maxiter; iter++) {

                // Initialize particles (NOT timed)
                initializepoints(&points);

                // Reset mesh (NOT timed)
                for (int i = 0; i < GRID_X * GRID_Y; i++)
                    mesh_value[i] = 0.0;

                // Time ONLY interpolation
                double start = omp_get_wtime();
                interpolation(mesh_value, &points);
                double end = omp_get_wtime();

                total_interp_time += (end - start);
            }

            // Write to CSV
            fprintf(fp, "%d,%d,%d,%lld,%lf\n",
                    config+1, NX, NY,
                    (long long)NUM_Points,
                    total_interp_time);

            printf("Total Interp Time = %lf seconds\n",
                   total_interp_time);

            free(mesh_value);
            free(points.x);
            free(points.y);
        }
    }

    fclose(fp);

    printf("\nAll Experiment 01 cases completed.\n");
    printf("Results stored in exp1_results.csv\n");

    return 0;
}
