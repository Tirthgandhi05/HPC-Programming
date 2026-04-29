#ifndef INIT_H
#define INIT_H

#include <stdio.h>

/* -------------------------------------------------------
 * Point structure
 * active: 1 = in domain,  0 = left domain (inactive)
 * ------------------------------------------------------- */
typedef struct {
    double x, y;
    int    active;
} Points;

/* -------------------------------------------------------
 * Global simulation parameters (defined in main.cpp)
 * ------------------------------------------------------- */
extern int    GRID_X, GRID_Y, NX, NY;
extern int    NUM_Points, Maxiter;
extern double dx, dy;

/* -------------------------------------------------------
 * I/O and initialisation
 * ------------------------------------------------------- */
void read_points(FILE *file, Points *points, int n);
void initializepoints(Points *points, int n);

#endif /* INIT_H */
