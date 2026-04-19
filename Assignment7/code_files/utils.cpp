#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cfloat>
#include <omp.h>

#include "utils.h"

double min_val, max_val;

static double *private_mesh_block = nullptr;
static int private_mesh_nthreads = 0;
static int private_mesh_size = 0;

void alloc_private_meshes(int nthreads) {
    const int mesh_size = GRID_X * GRID_Y;

    if (private_mesh_block != nullptr &&
        private_mesh_nthreads == nthreads &&
        private_mesh_size == mesh_size) {
        return;
    }

    free_private_meshes();

    private_mesh_block = static_cast<double *>(
        std::calloc(static_cast<std::size_t>(nthreads) * static_cast<std::size_t>(mesh_size),
                    sizeof(double))
    );

    if (!private_mesh_block) {
        std::fprintf(stderr, "Error: failed to allocate private meshes.\n");
        std::exit(EXIT_FAILURE);
    }

    private_mesh_nthreads = nthreads;
    private_mesh_size = mesh_size;
}

void free_private_meshes() {
    std::free(private_mesh_block);
    private_mesh_block = nullptr;
    private_mesh_nthreads = 0;
    private_mesh_size = 0;
}

static inline void scatter_one_point(double *priv_mesh, const Points &pt) {
    const double px = pt.x;
    const double py = pt.y;

    int col = static_cast<int>(px / dx);
    int row = static_cast<int>(py / dy);

    if (col >= NX) col = NX - 1;
    if (row >= NY) row = NY - 1;

    const double lx = px - static_cast<double>(col) * dx;
    const double ly = py - static_cast<double>(row) * dy;

    const double w_bl = (dx - lx) * (dy - ly);
    const double w_br = lx * (dy - ly);
    const double w_tl = (dx - lx) * ly;
    const double w_tr = lx * ly;

    const int idx_bl = row * GRID_X + col;
    const int idx_br = idx_bl + 1;
    const int idx_tl = idx_bl + GRID_X;
    const int idx_tr = idx_tl + 1;

    priv_mesh[idx_bl] += w_bl;
    priv_mesh[idx_br] += w_br;
    priv_mesh[idx_tl] += w_tl;
    priv_mesh[idx_tr] += w_tr;
}

void interpolation(double *mesh_value, Points *points) {
    const int mesh_size = GRID_X * GRID_Y;
    const int nthreads = omp_get_max_threads();

    alloc_private_meshes(nthreads);

    std::memset(mesh_value, 0, static_cast<std::size_t>(mesh_size) * sizeof(double));

    #pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        double *priv_mesh = private_mesh_block + static_cast<std::size_t>(tid) * mesh_size;

        std::memset(priv_mesh, 0, static_cast<std::size_t>(mesh_size) * sizeof(double));

        #pragma omp for schedule(static)
        for (int p = 0; p < NUM_Points; ++p) {
            if (!points[p].is_void) {
                scatter_one_point(priv_mesh, points[p]);
            }
        }

        #pragma omp for schedule(static)
        for (int i = 0; i < mesh_size; ++i) {
            double sum = 0.0;
            for (int t = 0; t < nthreads; ++t) {
                sum += private_mesh_block[static_cast<std::size_t>(t) * mesh_size + i];
            }
            mesh_value[i] = sum;
        }
    }
}

void normalization(double *mesh_value) {
    const int mesh_size = GRID_X * GRID_Y;
    double g_min = DBL_MAX;
    double g_max = -DBL_MAX;

    #pragma omp parallel for schedule(static) reduction(min:g_min) reduction(max:g_max)
    for (int i = 0; i < mesh_size; ++i) {
        const double v = mesh_value[i];
        if (v < g_min) g_min = v;
        if (v > g_max) g_max = v;
    }

    min_val = g_min;
    max_val = g_max;

    const double range = max_val - min_val;
    if (range <= 1e-15) {
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < mesh_size; ++i) {
            mesh_value[i] = 0.0;
        }
        return;
    }

    const double scale = 2.0 / range;
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < mesh_size; ++i) {
        mesh_value[i] = (mesh_value[i] - min_val) * scale - 1.0;
    }
}

static inline double interpolate_from_mesh(const double *mesh_value, double px, double py) {
    int col = static_cast<int>(px / dx);
    int row = static_cast<int>(py / dy);

    if (col >= NX) col = NX - 1;
    if (row >= NY) row = NY - 1;

    const double lx = px - static_cast<double>(col) * dx;
    const double ly = py - static_cast<double>(row) * dy;

    const double w_bl = (dx - lx) * (dy - ly);
    const double w_br = lx * (dy - ly);
    const double w_tl = (dx - lx) * ly;
    const double w_tr = lx * ly;

    const int idx_bl = row * GRID_X + col;
    const int idx_br = idx_bl + 1;
    const int idx_tl = idx_bl + GRID_X;
    const int idx_tr = idx_tl + 1;

    return w_bl * mesh_value[idx_bl] +
           w_br * mesh_value[idx_br] +
           w_tl * mesh_value[idx_tl] +
           w_tr * mesh_value[idx_tr];
}

void mover(double *mesh_value, Points *points) {
    #pragma omp parallel for schedule(static)
    for (int p = 0; p < NUM_Points; ++p) {
        if (points[p].is_void) continue;

        const double px = points[p].x;
        const double py = points[p].y;

        const double Fi = interpolate_from_mesh(mesh_value, px, py);

        const double x_new = px + Fi * dx;
        const double y_new = py + Fi * dy;

        if (x_new < 0.0 || x_new > 1.0 || y_new < 0.0 || y_new > 1.0) {
            points[p].is_void = true;
        } else {
            points[p].x = x_new;
            points[p].y = y_new;
        }
    }
}

void denormalization(double *mesh_value) {
    const int mesh_size = GRID_X * GRID_Y;
    const double range = max_val - min_val;

    if (range <= 1e-15) {
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < mesh_size; ++i) {
            mesh_value[i] = min_val;
        }
        return;
    }

    const double half_range = 0.5 * range;
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < mesh_size; ++i) {
        mesh_value[i] = (mesh_value[i] + 1.0) * half_range + min_val;
    }
}

long long int void_count(Points *points) {
    long long int voids = 0;
    #pragma omp parallel for schedule(static) reduction(+:voids)
    for (int i = 0; i < NUM_Points; ++i) {
        voids += static_cast<long long int>(points[i].is_void);
    }
    return voids;
}

void save_mesh(double *mesh_value) {
    FILE *fd = std::fopen("Mesh.out", "w");
    if (!fd) {
        std::fprintf(stderr, "Error creating Mesh.out\n");
        std::exit(EXIT_FAILURE);
    }

    for (int i = 0; i < GRID_Y; ++i) {
        for (int j = 0; j < GRID_X; ++j) {
            std::fprintf(fd, "%.6f ", mesh_value[static_cast<std::size_t>(i) * GRID_X + j]);
        }
        std::fprintf(fd, "\n");
    }

    std::fclose(fd);
}
