/*=============================================================
 * utils.cpp  –  Serial implementations of mesh operations
 *
 * interpolation   : particles -> mesh  (CIC bilinear scatter)
 * normalise_mesh  : mesh values -> [-1,1]
 * mover           : mesh -> particles  (reverse CIC gather)
 * denormalise_mesh: undo normalisation
 * save_mesh       : write mesh to Mesh.out
 *=============================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "utils.h"
#include "init.h"

/* -----------------------------------------------------------
 * interpolation  –  Cloud-in-Cell (bilinear) scatter
 *
 * For each ACTIVE particle (fi = 1):
 *   gx = floor(px/dx),   gy = floor(py/dy)
 *   lx = px - gx*dx,     ly = py - gy*dy
 *   Weights:
 *     a1 = (dx-lx)*(dy-ly)  -> node (gx,   gy  )
 *     a2 = lx*(dy-ly)        -> node (gx+1, gy  )
 *     a3 = (dx-lx)*ly        -> node (gx,   gy+1)
 *     a4 = lx*ly             -> node (gx+1, gy+1)
 *
 * Caller must zero mesh_value before calling.
 * ----------------------------------------------------------- */
void interpolation(double *mesh_value, Points *points, int n_pts)
{
    int p;
    for (p = 0; p < n_pts; p++) {

        if (!points[p].active) continue;

        double px = points[p].x;
        double py = points[p].y;

        int gx = (int)(px / dx);
        int gy = (int)(py / dy);

        /* Clamp to valid cell range */
        if (gx >= NX) gx = NX - 1;
        if (gy >= NY) gy = NY - 1;
        if (gx < 0)   gx = 0;
        if (gy < 0)   gy = 0;

        double lx = px - gx * dx;
        double ly = py - gy * dy;

        double a1 = (dx - lx) * (dy - ly);
        double a2 = lx        * (dy - ly);
        double a3 = (dx - lx) * ly;
        double a4 = lx        * ly;

        int base = gy * GRID_X + gx;
        mesh_value[base]              += a1;
        mesh_value[base + 1]          += a2;
        mesh_value[base + GRID_X]     += a3;
        mesh_value[base + GRID_X + 1] += a4;
    }
}

/* -----------------------------------------------------------
 * normalise_mesh  –  linearly map mesh values to [-1, 1]
 *   F_norm = 2*(F - Fmin)/(Fmax-Fmin) - 1
 * ----------------------------------------------------------- */
void normalise_mesh(double *mesh_value, double *out_min, double *out_max)
{
    int n = GRID_X * GRID_Y;
    double fmin =  DBL_MAX;
    double fmax = -DBL_MAX;
    int i;

    for (i = 0; i < n; i++) {
        if (mesh_value[i] < fmin) fmin = mesh_value[i];
        if (mesh_value[i] > fmax) fmax = mesh_value[i];
    }

    *out_min = fmin;
    *out_max = fmax;

    double range = fmax - fmin;
    if (range < 1e-15) {
        memset(mesh_value, 0, n * sizeof(double));
        return;
    }

    double inv_range = 2.0 / range;
    for (i = 0; i < n; i++)
        mesh_value[i] = (mesh_value[i] - fmin) * inv_range - 1.0;
}

/* -----------------------------------------------------------
 * mover  –  reverse CIC gather + position update
 *
 * For each ACTIVE particle:
 *   F = sum_k(weight_k * mesh[k])   (same weights as scatter)
 *   x_new = x + F*dx
 *   y_new = y + F*dy
 *   If outside [0,1]^2 -> mark inactive
 * ----------------------------------------------------------- */
void mover(double *mesh_value, Points *points, int n_pts)
{
    int p;
    for (p = 0; p < n_pts; p++) {

        if (!points[p].active) continue;

        double px = points[p].x;
        double py = points[p].y;

        int gx = (int)(px / dx);
        int gy = (int)(py / dy);

        if (gx >= NX) gx = NX - 1;
        if (gy >= NY) gy = NY - 1;
        if (gx < 0)   gx = 0;
        if (gy < 0)   gy = 0;

        double lx = px - gx * dx;
        double ly = py - gy * dy;

        double a1 = (dx - lx) * (dy - ly);
        double a2 = lx        * (dy - ly);
        double a3 = (dx - lx) * ly;
        double a4 = lx        * ly;

        int base = gy * GRID_X + gx;
        double F = a1 * mesh_value[base]
                 + a2 * mesh_value[base + 1]
                 + a3 * mesh_value[base + GRID_X]
                 + a4 * mesh_value[base + GRID_X + 1];

        double x_new = px + F * dx;
        double y_new = py + F * dy;

        if (x_new < 0.0 || x_new > 1.0 ||
            y_new < 0.0 || y_new > 1.0) {
            points[p].active = 0;
        } else {
            points[p].x = x_new;
            points[p].y = y_new;
        }
    }
}

/* -----------------------------------------------------------
 * denormalise_mesh  –  undo normalise_mesh
 *   F = (F_norm + 1) * (Fmax-Fmin)/2 + Fmin
 * ----------------------------------------------------------- */
void denormalise_mesh(double *mesh_value, double mesh_min, double mesh_max)
{
    int n = GRID_X * GRID_Y;
    double range = mesh_max - mesh_min;
    int i;

    if (range < 1e-15) {
        for (i = 0; i < n; i++) mesh_value[i] = mesh_min;
        return;
    }

    double half_range = range * 0.5;
    for (i = 0; i < n; i++)
        mesh_value[i] = (mesh_value[i] + 1.0) * half_range + mesh_min;
}

/* -----------------------------------------------------------
 * save_mesh  –  write GRID_Y rows x GRID_X cols to Mesh.out
 * ----------------------------------------------------------- */
void save_mesh(double *mesh_value)
{
    FILE *fd = fopen("Mesh.out", "w");
    if (!fd) {
        printf("Error: cannot create Mesh.out\n");
        return;
    }

    int i, j;
    for (i = 0; i < GRID_Y; i++) {
        for (j = 0; j < GRID_X - 1; j++)
            fprintf(fd, "%.6lf ", mesh_value[i * GRID_X + j]);
        fprintf(fd, "%.6lf\n", mesh_value[i * GRID_X + GRID_X - 1]);
    }
    fclose(fd);
}
