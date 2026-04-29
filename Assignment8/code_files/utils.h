#ifndef UTILS_H
#define UTILS_H

#include "init.h"

/* -------------------------------------------------------
 * Forward interpolation  (particles → mesh)
 *
 * Uses bilinear (Cloud-in-Cell) weighting.
 * mesh_value MUST be zeroed by the caller before calling.
 * n_pts  = number of particles to process (allows partial
 *          arrays for MPI local subsets).
 * ------------------------------------------------------- */
void interpolation(double *mesh_value, Points *points, int n_pts);

/* -------------------------------------------------------
 * Normalise mesh values linearly to [-1, 1].
 * Stores original min/max so denormalise() can undo it.
 * ------------------------------------------------------- */
void normalise_mesh(double *mesh_value, double *out_min, double *out_max);

/* -------------------------------------------------------
 * Reverse interpolation  (mesh → particles)  "mover"
 *
 * Reads the normalised grid field back to each particle,
 * updates positions:
 *     x_new = x + F * dx
 *     y_new = y + F * dy
 * Particles that leave [0,1]x[0,1] are marked inactive.
 * ------------------------------------------------------- */
void mover(double *mesh_value, Points *points, int n_pts);

/* -------------------------------------------------------
 * Denormalise mesh back to original value range.
 * ------------------------------------------------------- */
void denormalise_mesh(double *mesh_value, double mesh_min, double mesh_max);

/* -------------------------------------------------------
 * Write mesh to "Mesh.out"
 * ------------------------------------------------------- */
void save_mesh(double *mesh_value);

#endif /* UTILS_H */
