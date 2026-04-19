#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <string>
#include <omp.h>

#include "init.h"
#include "utils.h"

// Global variable definitions
int GRID_X, GRID_Y, NX, NY;
int NUM_Points, Maxiter;
double dx, dy;

struct RunStats {
    double interp_time = 0.0;
    double norm_time = 0.0;
    double move_time = 0.0;
    double denorm_time = 0.0;
    double total_time = 0.0;
    long long voids = 0;
};

static void read_header_and_points(const char *input_file,
                                  std::vector<Points> &points_out)
{
    FILE *file = std::fopen(input_file, "rb");
    if (!file) {
        std::fprintf(stderr, "Error opening input file: %s\n", input_file);
        std::exit(EXIT_FAILURE);
    }

    if (std::fread(&NX, sizeof(int), 1, file) != 1 ||
        std::fread(&NY, sizeof(int), 1, file) != 1 ||
        std::fread(&NUM_Points, sizeof(int), 1, file) != 1 ||
        std::fread(&Maxiter, sizeof(int), 1, file) != 1) {
        std::fprintf(stderr, "Error: input file header is incomplete.\n");
        std::exit(EXIT_FAILURE);
    }

    GRID_X = NX + 1;
    GRID_Y = NY + 1;
    dx = 1.0 / static_cast<double>(NX);
    dy = 1.0 / static_cast<double>(NY);

    points_out.resize(NUM_Points);
    read_points(file, points_out.data());

    std::fclose(file);
}

static RunStats run_simulation(int nthreads,
                               const std::vector<Points> &initial_points,
                               std::vector<Points> &work_points,
                               std::vector<double> &mesh_value,
                               bool save_output_mesh)
{
    omp_set_dynamic(0);
    omp_set_num_threads(nthreads);

    work_points = initial_points;
    std::fill(mesh_value.begin(), mesh_value.end(), 0.0);

    RunStats stats;

    for (int iter = 0; iter < Maxiter; ++iter) {
        const double t0 = omp_get_wtime();
        interpolation(mesh_value.data(), work_points.data());
        const double t1 = omp_get_wtime();

        normalization(mesh_value.data());
        const double t2 = omp_get_wtime();

        mover(mesh_value.data(), work_points.data());
        const double t3 = omp_get_wtime();

        denormalization(mesh_value.data());
        const double t4 = omp_get_wtime();

        stats.interp_time += (t1 - t0);
        stats.norm_time   += (t2 - t1);
        stats.move_time   += (t3 - t2);
        stats.denorm_time += (t4 - t3);
    }

    stats.total_time = stats.interp_time + stats.norm_time + stats.move_time + stats.denorm_time;
    stats.voids = void_count(work_points.data());

    if (save_output_mesh) {
        save_mesh(mesh_value.data());
    }

    return stats;
}

static void print_stats(const char *label, int nthreads, const RunStats &stats,
                        double baseline_time = 0.0)
{
    std::printf("\n--- %s (%d thread%s) ---\n", label, nthreads, (nthreads == 1 ? "" : "s"));
    std::printf("Total Interpolation Time  = %.6f seconds\n", stats.interp_time);
    std::printf("Total Normalization Time  = %.6f seconds\n", stats.norm_time);
    std::printf("Total Mover Time          = %.6f seconds\n", stats.move_time);
    std::printf("Total Denormalization Time= %.6f seconds\n", stats.denorm_time);
    std::printf("Total Algorithm Time      = %.6f seconds\n", stats.total_time);
    std::printf("Total Number of Voids     = %lld\n", stats.voids);

    if (baseline_time > 0.0) {
        const double speedup = baseline_time / stats.total_time;
        const double eff = speedup / static_cast<double>(nthreads);
        std::printf("Speedup                   = %.4fx\n", speedup);
        std::printf("Parallel Efficiency       = %.4f\n", eff);
    }
}

static void benchmark_all(const std::vector<Points> &initial_points,
                          std::vector<Points> &work_points,
                          std::vector<double> &mesh_value)
{
    const int thread_counts[] = {1, 2, 4, 8, 16};
    const int num_thread_counts = 5;

    FILE *csv = std::fopen("results.csv", "a");
    if (!csv) {
        std::fprintf(stderr, "Error opening results.csv for writing.\n");
        std::exit(EXIT_FAILURE);
    }

    std::fprintf(csv, "Config,NX,NY,Points,Maxiter,Threads,TotalTime,InterpTime,NormTime,MoverTime,DenormTime,Speedup,Efficiency,Voids\n");

    RunStats baseline = run_simulation(1, initial_points, work_points, mesh_value, false);
    print_stats("Serial baseline", 1, baseline);

    std::fprintf(csv, "NX%d_NY%d_P%d,%d,%d,%d,%d,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.0f\n",
                 NX, NY, NUM_Points, NX, NY, NUM_Points, Maxiter, 1,
                 baseline.total_time, baseline.interp_time, baseline.norm_time,
                 baseline.move_time, baseline.denorm_time, 1.0, 1.0,
                 static_cast<double>(baseline.voids));

    for (int i = 1; i < num_thread_counts; ++i) {
        const int t = thread_counts[i];
        RunStats stats = run_simulation(t, initial_points, work_points, mesh_value, false);
        print_stats("Parallel", t, stats, baseline.total_time);

        const double speedup = baseline.total_time / stats.total_time;
        const double efficiency = speedup / static_cast<double>(t);

        std::fprintf(csv, "NX%d_NY%d_P%d,%d,%d,%d,%d,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.0f\n",
                     NX, NY, NUM_Points, NX, NY, NUM_Points, Maxiter, t,
                     stats.total_time, stats.interp_time, stats.norm_time,
                     stats.move_time, stats.denorm_time, speedup, efficiency,
                     static_cast<double>(stats.voids));
    }

    std::fclose(csv);

    // Keep the last computed mesh on disk so the user can inspect it.
    save_mesh(mesh_value.data());

    std::printf("\nBenchmark results written to results.csv\n");
}

int main(int argc, char **argv) {
    if (argc < 2) {
        std::printf("Usage: %s <input_file> [--bench]\n", argv[0]);
        std::printf("Run normally with OMP_NUM_THREADS set to 1,2,4,8,16.\n");
        return 1;
    }

    bool benchmark = false;
    for (int i = 2; i < argc; ++i) {
        if (std::strcmp(argv[i], "--bench") == 0) {
            benchmark = true;
        }
    }

    std::vector<Points> initial_points;
    read_header_and_points(argv[1], initial_points);

    std::vector<Points> work_points(static_cast<std::size_t>(NUM_Points));
    std::vector<double> mesh_value(static_cast<std::size_t>(GRID_X) * GRID_Y, 0.0);

    std::printf("Grid: %d x %d | Particles: %d | Iterations: %d\n",
                NX, NY, NUM_Points, Maxiter);

    if (benchmark) {
        benchmark_all(initial_points, work_points, mesh_value);
        free_private_meshes();
        return 0;
    }

    int nthreads = omp_get_max_threads();
    if (nthreads < 1) nthreads = 1;

    RunStats stats = run_simulation(nthreads, initial_points, work_points, mesh_value, true);
    print_stats("Run complete", nthreads, stats);

    std::printf("\nOutput written to Mesh.out\n");
    free_private_meshes();
    return 0;
}
