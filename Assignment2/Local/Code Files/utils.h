#ifndef UTILS_H
#define UTILS_H
#include <time.h>

void matrix_multiplication(double** A, double** B, double** C, int N);
void transpose(double** B, double** BT, int N);
void transposed_matrix_multiplication(double** A, double** B, double** C, int N);
void block_matrix_multiplication(double** A, double** B, double** C, int N);

#endif
