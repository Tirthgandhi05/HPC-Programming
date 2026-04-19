#ifndef UTILS_H
#define UTILS_H
#include "init.h"

void interpolation(double *mesh_value, Points *points);
void parallel_interpolation(double *mesh_value, Points *points, double *all_local_meshes, int num_threads);
void save_mesh(double *mesh_value);

#endif
