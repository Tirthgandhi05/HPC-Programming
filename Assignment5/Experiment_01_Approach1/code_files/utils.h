#ifndef UTILS_H
#define UTILS_H

#include "init.h"

void interpolation(double *mesh_value, Points *points);

void mover_immediate(Points *points, double deltaX, double deltaY);
void mover_deferred(Points *points, double deltaX, double deltaY);

void save_mesh(double *mesh_value);

#endif