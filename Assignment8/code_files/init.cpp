#include <stdio.h>
#include <stdlib.h>
#include "init.h"

void read_points(FILE *file, Points *points, int n)
{
    int i;
    for (i = 0; i < n; i++) {
        fread(&points[i].x, sizeof(double), 1, file);
        fread(&points[i].y, sizeof(double), 1, file);
        points[i].active = 1;
    }
}

void initializepoints(Points *points, int n)
{
    int i;
    for (i = 0; i < n; i++) {
        points[i].x      = (double)rand() / RAND_MAX;
        points[i].y      = (double)rand() / RAND_MAX;
        points[i].active = 1;
    }
}
