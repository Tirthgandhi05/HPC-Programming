#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "utils.h"
#include "init.h"

int main()
{

    int grids[3][2] = {{250, 100}, {500, 200}, {1000, 400}};
    int threads_list[4] = {2, 4, 8, 16};

    int N = 14000000;
    int Maxiter = 10;

    for (int g = 0; g < 3; g++)
    {

        int NX = grids[g][0];
        int NY = grids[g][1];
        int GRID_X = NX + 1;
        double dx = 1.0 / NX;
        double dy = 1.0 / NY;

        printf("=========== GRID %dx%d ===========\n", NX, NY);

        Points p_im, p_def;
        p_im.x = (double *)malloc(N * sizeof(double));
        p_im.y = (double *)malloc(N * sizeof(double));
        p_def.x = (double *)malloc(N * sizeof(double));
        p_def.y = (double *)malloc(N * sizeof(double));

        // ✅ CORRECT SERIAL BASELINE
        initializepoints(&p_im, N);

        double serial_im = 0;
        for (int i = 0; i < Maxiter; i++)
        {
            double t = omp_get_wtime();
            mover_immediate_serial(&p_im, N, dx, dy);
            serial_im += (omp_get_wtime() - t);
        }
        serial_im /= Maxiter;

        initializepoints(&p_def, N);

        double serial_def = 0;
        for (int i = 0; i < Maxiter; i++)
        {
            double t = omp_get_wtime();
            mover_deferred_serial(&p_def, N, dx, dy);
            serial_def += (omp_get_wtime() - t);
        }
        serial_def /= Maxiter;

        for (int t = 0; t < 4; t++)
        {

            int threads = threads_list[t];

            // ✅ FAIR DATA RESET
            initializepoints(&p_im, N);
            for (int i = 0; i < N; i++)
            {
                p_def.x[i] = p_im.x[i];
                p_def.y[i] = p_im.y[i];
            }

            double mover_im = 0, mover_def = 0;

            for (int iter = 0; iter < Maxiter; iter++)
            {

                double t1 = omp_get_wtime();
                mover_immediate_parallel(&p_im, N, dx, dy, threads);
                mover_im += (omp_get_wtime() - t1);

                double t2 = omp_get_wtime();
                mover_deferred_parallel(&p_def, N, dx, dy, threads);
                mover_def += (omp_get_wtime() - t2);
            }

            mover_im /= Maxiter;
            mover_def /= Maxiter;

            double sp_im = serial_im / mover_im;
            double sp_def = serial_def / mover_def;

            printf("[Threads %d] Immediate: %.4fs (%.2fx) | Deferred: %.4fs (%.2fx)\n",
                   threads, mover_im, sp_im, mover_def, sp_def);
        }

        free(p_im.x);
        free(p_im.y);
        free(p_def.x);
        free(p_def.y);

        printf("\n");
    }

    return 0;
}
