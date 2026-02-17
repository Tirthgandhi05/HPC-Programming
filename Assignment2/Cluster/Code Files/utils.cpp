#include <math.h>
#include "utils.h"
#include <cstdlib>

// Problem 01
void matrix_multiplication(double** A, double** B, double** C, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            C[i][j] = 0.0;
            for (int k = 0; k < N; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

void transpose(double** B, double** BT, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            BT[j][i] = B[i][j];
        }
    }
}

void transposed_matrix_multiplication(double** A, double** B, double** C, int N) {
    double** BT = (double**)malloc(N * sizeof(double*));
    for (int i = 0; i < N; i++)
        BT[i] = (double*)malloc(N * sizeof(double));

    transpose(B, BT, N);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            C[i][j] = 0.0;
            for (int k = 0; k < N; k++) {
                C[i][j] += A[i][k] * BT[j][k];
            }
        }
    }

    for (int i = 0; i < N; i++)
        free(BT[i]);
    free(BT);
}

void block_matrix_multiplication(double** A, double** B, double** C, int N) {

    // Choose block size (cache-friendly)
    int Bsize;
    if (N >= 64)
        Bsize = 32;
    else
        Bsize = N;

    // Initialize result matrix
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            C[i][j] = 0.0;

    // Blocked multiplication
    for (int ii = 0; ii < N; ii += Bsize) {
        for (int jj = 0; jj < N; jj += Bsize) {
            for (int kk = 0; kk < N; kk += Bsize) {

                for (int i = ii; i < ii + Bsize && i < N; i++) {
                    for (int j = jj; j < jj + Bsize && j < N; j++) {
                        double sum = C[i][j];
                        for (int k = kk; k < kk + Bsize && k < N; k++) {
                            sum += A[i][k] * B[k][j];
                        }
                        C[i][j] = sum;
                    }
                }

            }
        }
    }
}

