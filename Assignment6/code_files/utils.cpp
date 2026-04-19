#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include "utils.h"
#include "init.h"

void interpolation(double *mesh_value, Points *points) {
    double inv_dx = 1.0 / dx;
    double inv_dy = 1.0 / dy;

    for (int p = 0; p < NUM_Points; p++) {
        double x = points[p].x;
        double y = points[p].y;

        int i = (int)(x * inv_dx);
        int j = (int)(y * inv_dy);

        if (i == NX) i = NX - 1;
        if (j == NY) j = NY - 1;

        double dx_local = x - i * dx;
        double dy_local = y - j * dy;

        int base = j * GRID_X + i;

        mesh_value[base]              += (dx - dx_local) * (dy - dy_local);
        mesh_value[base + 1]          += dx_local * (dy - dy_local);
        mesh_value[base + GRID_X]     += (dx - dx_local) * dy_local;
        mesh_value[base + GRID_X + 1] += dx_local * dy_local ;
    }
}

void parallel_interpolation(double *mesh_value, Points *points,
                            double *all_local_meshes, int num_threads) {

    double inv_dx = 1.0 / dx;
    double inv_dy = 1.0 / dy;

    int total_grid_points = GRID_X * GRID_Y;

    #pragma omp parallel num_threads(num_threads)
    {
        int tid = omp_get_thread_num();
        double *local_mesh = &all_local_meshes[tid * total_grid_points];

        // ✅ Reset local mesh
        memset(local_mesh, 0, total_grid_points * sizeof(double));

        // ✅ Parallel interpolation (no race)
        #pragma omp for schedule(static)
        for (int p = 0; p < NUM_Points; p++) {
            double x = points[p].x;
            double y = points[p].y;

            int i = (int)(x * inv_dx);
            int j = (int)(y * inv_dy);

            if (i == NX) i = NX - 1;
            if (j == NY) j = NY - 1;

            double dx_local = x - i * dx;
            double dy_local = y - j * dy;

            int base = j * GRID_X + i;

            local_mesh[base]              += (dx - dx_local) * (dy - dy_local);
            local_mesh[base + 1]          += dx_local * (dy - dy_local);
            local_mesh[base + GRID_X]     += (dx - dx_local) * dy_local;
            local_mesh[base + GRID_X + 1] += dx_local * dy_local;
        }

        // ✅ Reduction
        #pragma omp for schedule(static)
        for (int k = 0; k < total_grid_points; k++) {
            double sum = 0.0;
            for (int t = 0; t < num_threads; t++) {
                sum += all_local_meshes[t * total_grid_points + k];
            }
            mesh_value[k] = sum;
        }
    }
}

void save_mesh(double *mesh_value) {
    FILE *fd = fopen("Mesh.out", "w");
    if (!fd) exit(1);

    for (int i = 0; i < GRID_Y; i++) {
        for (int j = 0; j < GRID_X; j++) {
            fprintf(fd, "%.6lf ", mesh_value[i * GRID_X + j]);
        }
        fprintf(fd, "\n");
    }

    fclose(fd);
}
