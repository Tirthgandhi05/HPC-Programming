#ifndef UTILS_H
#define UTILS_H

void dummy(int x);

void kernel_copy(double *x, double *y, int Np);
void kernel_scale(double *x, double *y, double a, int Np);
void kernel_add(double *x, double *y, double *S, int Np);
void kernel_triad(double *x, double *y, double *S, double a, int Np);

#endif

