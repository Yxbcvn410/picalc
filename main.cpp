#include <iostream>
#include <fstream>
#include <cmath>

#include <mpi.h>
// Электронная почта: carpson@mail.ru

using namespace std;

double pi_approx_simple(long N, int k, int rk) {
    double answer = (rk == 0 ? -0.5 : 0);
    double nn = 1 / (double) N / (double) N;
    for (long i = rk; i < N; i += k) {
        answer += sqrt(1 - (double) (i * i) * nn);
    }
    return answer * 4 / (double) N;
}

int main(int argc, char **argv) {
    int r = MPI_Init(&argc, &argv);
    if (r != 0) {
        printf("MPI init failed\n");
        return r;
    }

    int k, rank;
    long N;
    MPI_Comm_size(MPI_COMM_WORLD, &k);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        ifstream fs("N.dat", ios::in);
        fs >> N;
    }

    MPI_Bcast(&N, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    double result = pi_approx_simple(N, k, rank);
    double f = 0;

    MPI_Reduce(&result, &f, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        printf("%.17f\n", f);
    }

    MPI_Finalize();
}
