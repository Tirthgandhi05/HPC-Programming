/*
 * MPI_hybrid.c  –  Hybrid MPI+OpenMP PIC Interpolation + Mover
 *
 * Correctness goals:
 *   - particle decomposition across MPI ranks
 *   - thread-private accumulation for interpolation (no races)
 *   - MPI_Allreduce for the global mesh
 *   - mover on each rank's own particles
 *   - fair timing using the slowest rank (MPI_MAX)
 *
 * Build:
 *   mpicc -O3 -fopenmp -march=native MPI_hybrid_corrected.c -o mpi_hybrid -lm
 *
 * Run:
 *   mpirun -np <ranks> --hostfile sources.txt --map-by node \
 *          ./mpi_hybrid <input.bin> <omp_threads_per_rank>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <omp.h>
#include <mpi.h>

typedef struct { double x, y; int active; int _pad; } Points;

static int    GRID_X, GRID_Y, NX, NY, NUM_Points, NTHR_C, Maxiter;
static double dx, dy;
static MPI_Datatype MPI_POINTS_T;

static void create_mpi_points_type(void)
{
    Points dummy;
    MPI_Aint base, d[3];
    MPI_Get_address(&dummy,        &base);
    MPI_Get_address(&dummy.x,      &d[0]);
    MPI_Get_address(&dummy.y,      &d[1]);
    MPI_Get_address(&dummy.active, &d[2]);
    d[0] -= base; d[1] -= base; d[2] -= base;

    int bl[3]          = {1, 1, 1};
    MPI_Datatype tp[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_INT};
    MPI_Type_create_struct(3, bl, d, tp, &MPI_POINTS_T);
    MPI_Type_commit(&MPI_POINTS_T);
}

static inline void cic(double px, double py,
                       int *p1, int *p2, int *p3, int *p4,
                       double *a1, double *a2, double *a3, double *a4)
{
    int gx = (int)(px / dx);
    int gy = (int)(py / dy);

    if (gx >= NX) gx = NX - 1;
    if (gy >= NY) gy = NY - 1;
    if (gx < 0)   gx = 0;
    if (gy < 0)   gy = 0;

    double lx = px - gx * dx;
    double ly = py - gy * dy;

    *p1 =  gy      * GRID_X + gx;
    *p2 =  gy      * GRID_X + (gx + 1);
    *p3 = (gy + 1) * GRID_X + gx;
    *p4 = (gy + 1) * GRID_X + (gx + 1);

    *a1 = (dx - lx) * (dy - ly);
    *a2 =  lx       * (dy - ly);
    *a3 = (dx - lx) *  ly;
    *a4 =  lx       *  ly;
}

static void interpolation(double * restrict lmesh,
                          double * restrict priv,
                          Points * restrict lpts,
                          int ln, int msz)
{
    memset(priv, 0, (size_t)NTHR_C * msz * sizeof(double));

    #pragma omp parallel num_threads(NTHR_C)
    {
        int tid = omp_get_thread_num();
        double *pm = priv + (size_t)tid * msz;

        #pragma omp for schedule(static)
        for (int i = 0; i < ln; i++) {
            if (!lpts[i].active) continue;

            int p1, p2, p3, p4;
            double a1, a2, a3, a4;
            cic(lpts[i].x, lpts[i].y, &p1, &p2, &p3, &p4, &a1, &a2, &a3, &a4);

            pm[p1] += a1;
            pm[p2] += a2;
            pm[p3] += a3;
            pm[p4] += a4;
        }

        #pragma omp for schedule(static)
        for (int j = 0; j < msz; j++) {
            double acc = 0.0;
            for (int t = 0; t < NTHR_C; t++) {
                acc += priv[(size_t)t * msz + j];
            }
            lmesh[j] = acc;
        }
    }
}

static void mover(double * restrict gmesh,
                  Points * restrict lpts,
                  int ln)
{
    #pragma omp parallel for num_threads(NTHR_C) schedule(static)
    for (int i = 0; i < ln; i++) {
        if (!lpts[i].active) continue;

        int p1, p2, p3, p4;
        double a1, a2, a3, a4;
        cic(lpts[i].x, lpts[i].y, &p1, &p2, &p3, &p4, &a1, &a2, &a3, &a4);

        double F = a1 * gmesh[p1] + a2 * gmesh[p2]
                 + a3 * gmesh[p3] + a4 * gmesh[p4];

        double nx = lpts[i].x + F * dx;
        double ny = lpts[i].y + F * dy;

        if (nx < 0.0 || nx > 1.0 || ny < 0.0 || ny > 1.0) {
            lpts[i].active = 0;
        } else {
            lpts[i].x = nx;
            lpts[i].y = ny;
        }
    }
}

static void readPoints(FILE *f, Points *pts, int n)
{
    for (int i = 0; i < n; i++) {
        fread(&pts[i].x, sizeof(double), 1, f);
        fread(&pts[i].y, sizeof(double), 1, f);
        pts[i].active = 1;
        pts[i]._pad   = 0;
    }
}

static void printmesh(double *m)
{
    FILE *fd = fopen("Mesh.out", "w");
    if (!fd) {
        fprintf(stderr, "Cannot create Mesh.out\n");
        exit(1);
    }

    for (int i = 0; i < GRID_Y; i++) {
        for (int j = 0; j < GRID_X - 1; j++) {
            fprintf(fd, "%lf ", m[i * GRID_X + j]);
        }
        fprintf(fd, "%lf\n", m[i * GRID_X + GRID_X - 1]);
    }
    fclose(fd);
}

int main(int argc, char **argv)
{
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc != 3) {
        if (rank == 0) {
            fprintf(stderr, "Usage: %s <input_file> <omp_threads>\n", argv[0]);
        }
        MPI_Finalize();
        return 1;
    }

    NTHR_C = atoi(argv[2]);
    if (NTHR_C < 1) NTHR_C = 1;

    create_mpi_points_type();

    if (rank == 0) {
        FILE *f = fopen(argv[1], "rb");
        if (!f) {
            fprintf(stderr, "Cannot open %s\n", argv[1]);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        fread(&NX,         sizeof(int), 1, f);
        fread(&NY,         sizeof(int), 1, f);
        fread(&NUM_Points,  sizeof(int), 1, f);
        fread(&Maxiter,     sizeof(int), 1, f);
        fclose(f);
    }

    int params[4] = {NX, NY, NUM_Points, Maxiter};
    MPI_Bcast(params, 4, MPI_INT, 0, MPI_COMM_WORLD);
    NX = params[0]; NY = params[1]; NUM_Points = params[2]; Maxiter = params[3];

    GRID_X = NX + 1;
    GRID_Y = NY + 1;
    dx = 1.0 / NX;
    dy = 1.0 / NY;
    int msz = GRID_X * GRID_Y;

    int base_c = NUM_Points / size;
    int rem    = NUM_Points % size;
    int *counts = (int *)malloc(size * sizeof(int));
    int *displs  = (int *)malloc(size * sizeof(int));
    if (!counts || !displs) {
        fprintf(stderr, "Rank %d: allocation failed\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    for (int i = 0; i < size; i++) {
        counts[i] = (i < rem) ? base_c + 1 : base_c;
        displs[i] = (i == 0) ? 0 : displs[i - 1] + counts[i - 1];
    }

    int ln = counts[rank];

    if (rank == 0) {
        printf("==========================================================\n");
        printf("  Hybrid MPI+OpenMP PIC\n");
        printf("==========================================================\n");
        printf("  NX=%d NY=%d NUM_Points=%d Maxiter=%d\n", NX, NY, NUM_Points, Maxiter);
        printf("  MPI_ranks=%d  OMP_threads=%d  Total_cores=%d\n", size, NTHR_C, size * NTHR_C);
    }

    Points *all_pts = NULL;
    if (rank == 0) {
        all_pts = (Points *)malloc(NUM_Points * sizeof(Points));
        if (!all_pts) {
            fprintf(stderr, "OOM on root\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        FILE *f = fopen(argv[1], "rb");
        if (!f) {
            fprintf(stderr, "Cannot open %s\n", argv[1]);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        int hdr[4];
        fread(hdr, sizeof(int), 4, f);
        readPoints(f, all_pts, NUM_Points);
        fclose(f);
    }

    Points *lpts = (Points *)malloc((size_t)ln * sizeof(Points));
    if (!lpts) {
        fprintf(stderr, "OOM lpts rank %d\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    MPI_Scatterv(all_pts, counts, displs, MPI_POINTS_T,
                 lpts, ln, MPI_POINTS_T, 0, MPI_COMM_WORLD);

    double *lmesh = (double *)calloc((size_t)msz, sizeof(double));
    double *gmesh = (double *)calloc((size_t)msz, sizeof(double));
    double *priv  = (double *)calloc((size_t)NTHR_C * msz, sizeof(double));
    if (!lmesh || !gmesh || !priv) {
        fprintf(stderr, "OOM mesh rank %d\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    double local_interp = 0.0;
    double local_mover  = 0.0;
    double local_total  = 0.0;
    double vmin_save = 0.0, vmax_save = 0.0;

    for (int iter = 0; iter < Maxiter; iter++) {
        MPI_Barrier(MPI_COMM_WORLD);
        double iter_start = MPI_Wtime();

        memset(lmesh, 0, (size_t)msz * sizeof(double));

        double t0 = MPI_Wtime();
        interpolation(lmesh, priv, lpts, ln, msz);
        double t1 = MPI_Wtime();
        local_interp += (t1 - t0);

        MPI_Allreduce(lmesh, gmesh, msz, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        double vmin = DBL_MAX, vmax = -DBL_MAX;
        #pragma omp parallel for num_threads(NTHR_C) reduction(min:vmin) reduction(max:vmax) schedule(static)
        for (int j = 0; j < msz; j++) {
            if (gmesh[j] < vmin) vmin = gmesh[j];
            if (gmesh[j] > vmax) vmax = gmesh[j];
        }

        vmin_save = vmin;
        vmax_save = vmax;

        double rng = vmax - vmin;
        if (rng < 1e-15) {
            memset(gmesh, 0, (size_t)msz * sizeof(double));
        } else {
            double inv = 2.0 / rng;
            #pragma omp parallel for num_threads(NTHR_C) schedule(static)
            for (int j = 0; j < msz; j++) {
                gmesh[j] = (gmesh[j] - vmin) * inv - 1.0;
            }
        }

        double t2 = MPI_Wtime();
        mover(gmesh, lpts, ln);
        double t3 = MPI_Wtime();
        local_mover += (t3 - t2);

        if (rng < 1e-15) {
            #pragma omp parallel for num_threads(NTHR_C) schedule(static)
            for (int j = 0; j < msz; j++) {
                gmesh[j] = vmin;
            }
        } else {
            double half = 0.5 * rng;
            double mid  = 0.5 * (vmax + vmin);
            #pragma omp parallel for num_threads(NTHR_C) schedule(static)
            for (int j = 0; j < msz; j++) {
                gmesh[j] = gmesh[j] * half + mid;
            }
        }

        local_total += (MPI_Wtime() - iter_start);
    }

    double global_interp = 0.0, global_mover = 0.0, global_total = 0.0;
    MPI_Reduce(&local_interp, &global_interp, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_mover,  &global_mover,  1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_total,  &global_total,  1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        printmesh(gmesh);
        printf("----------------------------------------------------------\n");
        printf("  Total interpolation time : %.6f s\n", global_interp);
        printf("  Total mover time         : %.6f s\n", global_mover);
        printf("  Total pipeline time      : %.6f s\n", global_total);
        printf("==========================================================\n");
        free(all_pts);
    }

    free(lpts);
    free(lmesh);
    free(gmesh);
    free(priv);
    free(counts);
    free(displs);
    MPI_Type_free(&MPI_POINTS_T);
    MPI_Finalize();
    return 0;
}
