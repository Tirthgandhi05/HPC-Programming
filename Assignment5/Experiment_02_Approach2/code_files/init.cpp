#include "utils.h"
#include <stdint.h>

void initializepoints(Points *points, int N) {

    uint32_t seed = 123456789;

    for (int i = 0; i < N; i++) {
        points->x[i] = fast_rand_unit(&seed);
        points->y[i] = fast_rand_unit(&seed);
    }
}
