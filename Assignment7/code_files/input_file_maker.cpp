#include <cstdio>
#include <cstdlib>
#include <ctime>

#define INPUT_FILENAME "input.bin"

static void generate_input_file(int NX, int NY, int NUM_Points, int Maxiter) {
    FILE *fp = std::fopen(INPUT_FILENAME, "wb");
    if (!fp) {
        std::perror("Error: unable to create input.bin");
        std::exit(EXIT_FAILURE);
    }

    std::fwrite(&NX, sizeof(int), 1, fp);
    std::fwrite(&NY, sizeof(int), 1, fp);
    std::fwrite(&NUM_Points, sizeof(int), 1, fp);
    std::fwrite(&Maxiter, sizeof(int), 1, fp);

    std::srand(static_cast<unsigned int>(std::time(nullptr)));

    for (int i = 0; i < NUM_Points; ++i) {
        double x = (double) std::rand() / (double) RAND_MAX;
        double y = (double) std::rand() / (double) RAND_MAX;
        std::fwrite(&x, sizeof(double), 1, fp);
        std::fwrite(&y, sizeof(double), 1, fp);
    }

    std::fclose(fp);
    std::printf("Input file '%s' generated successfully.\n", INPUT_FILENAME);
}

int main() {
    int NX, NY, NUM_Points, Maxiter;

    std::printf("Enter grid dimensions (NX NY): ");
    std::scanf("%d %d", &NX, &NY);

    std::printf("Enter number of particles: ");
    std::scanf("%d", &NUM_Points);

    std::printf("Enter number of iterations: ");
    std::scanf("%d", &Maxiter);

    generate_input_file(NX, NY, NUM_Points, Maxiter);
    return 0;
}
