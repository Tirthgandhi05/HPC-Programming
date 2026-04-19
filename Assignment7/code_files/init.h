#ifndef INIT_H
#define INIT_H

#include <cstdio>
#include <cstdbool>

struct Points {
    double x;
    double y;
    bool is_void;
};

extern int GRID_X, GRID_Y, NX, NY;
extern int NUM_Points, Maxiter;
extern double dx, dy;

void initializepoints(Points *points);
void read_points(FILE *file, Points *points);

#endif
