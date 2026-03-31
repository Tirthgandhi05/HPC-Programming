#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "init.h"

// Random particle initialization 
void initializepoints(Points *points) {
    for (int i = 0; i < NUM_Points; i++) {
        points->x[i] = (double) rand() / RAND_MAX;
        points->y[i] = (double) rand() / RAND_MAX;
    }
}