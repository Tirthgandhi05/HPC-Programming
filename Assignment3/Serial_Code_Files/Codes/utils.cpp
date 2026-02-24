#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include "init.h"
// Serial interpolation 
void interpolation(double *mesh_value, Points *points) {
	

    double inv_dx = 1.0 / dx;
    double inv_dy = 1.0 / dy;

    for (int p = 0; p < NUM_Points; p++) {
    

        double x = points[p].x;
        double y = points[p].y;

        // Convert to grid coordinates
        double gx = x * inv_dx;
        double gy = y * inv_dy;

        int i = (int)gx;
        int j = (int)gy;

        // Boundary safety
        if (i >= NX) i = NX - 1;
        if (j >= NY) j = NY - 1;

        // Local coordinates inside cell
        double dx_local = (x - i*dx);
        double dy_local = (y - j*dy);

        double one_minus_dx = dx - dx_local;
        double one_minus_dy = dy - dy_local;

        // Base index in 1D mesh array
        int base = j * GRID_X + i;

        mesh_value[base]                 += one_minus_dx * one_minus_dy;
        mesh_value[base + 1]             += dx_local * one_minus_dy;
        mesh_value[base + GRID_X]        += one_minus_dx * dy_local;
        mesh_value[base + GRID_X + 1]    += dx_local * dy_local;
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
