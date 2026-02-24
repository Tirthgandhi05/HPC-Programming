#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "init.h"

// Random particle initialization (optional)
void initializepoints(Points *points) {
    for (int i = 0; i < NUM_Points; i++) {
        points->x[i] = (double) rand() / RAND_MAX;
        points->y[i] = (double) rand() / RAND_MAX;
    }
}

// Read particle positions from binary file
void read_points(FILE *file, Points *points) {
    double temp[2];
    for (int i = 0; i < NUM_Points; i++) {
        fread(temp, sizeof(double), 2, file);
        points->x[i] = temp[0];
        points->y[i] = temp[1];
    }
}