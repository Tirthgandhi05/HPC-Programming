#include <cstdio>
#include <cstdlib>
#include <random>
#include "init.h"

void initializepoints(Points *points) {
    for (int i = 0; i < NUM_Points; i++) {
        points[i].x = (double) std::rand() / (double) RAND_MAX;
        points[i].y = (double) std::rand() / (double) RAND_MAX;
        points[i].is_void = false;
    }
}

void read_points(FILE *file, Points *points) {
    for (int i = 0; i < NUM_Points; i++) {
        if (std::fread(&points[i].x, sizeof(double), 1, file) != 1 ||
            std::fread(&points[i].y, sizeof(double), 1, file) != 1) {
            std::fprintf(stderr, "Error: unexpected end of input while reading points.\n");
            std::exit(EXIT_FAILURE);
        }
        points[i].is_void = false;
    }
}
