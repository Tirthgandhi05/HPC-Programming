#include <math.h>
#include "utils.h"

void kernel_copy(double *x, double *y, int Np) {
    for (int p = 0; p < Np; p++) {
        x[p] = y[p];

        if (((double)p) == 333.333)
            dummy(p);
    }
}

void kernel_scale(double *x, double *y, double a, int Np) {
    for (int p = 0; p < Np; p++) {
        x[p] = a * y[p];

        if (((double)p) == 333.333)
            dummy(p);
    }
}

void kernel_add(double *x, double *y, double *S, int Np) {
    for (int p = 0; p < Np; p++) {
        S[p] = x[p] + y[p];

        if (((double)p) == 333.333)
            dummy(p);
    }
}

void kernel_triad(double *x, double *y, double *S, double a, int Np) {
    for (int p = 0; p < Np; p++) {
        S[p] = x[p] + a * y[p];

        if (((double)p) == 333.333)
            dummy(p);
    }
}


void dummy(int x) {
    x = 10 * sin(x / 10.0);
}
