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
        for (int i = 1; i < k; ++i)
            MPI_Send(&N, sizeof(typeof(N)), MPI_LONG, i, 123, MPI_COMM_WORLD);
    } else {
        MPI_Recv(&N, sizeof(typeof(N)), MPI_LONG, 0, 123, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    double result = pi_approx_simple(N, k, rank);

    if (rank == 0) {
        double msg;
        for (int i = 1; i < k; ++i) {
            MPI_Recv(&msg, sizeof(typeof(msg)), MPI_DOUBLE, i, 234, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            result += msg;
        }
        printf("%.17f", result);
    } else {
        MPI_Send(&result, sizeof(typeof(result)), MPI_DOUBLE, 0, 234, MPI_COMM_WORLD);
    }

    MPI_Finalize();
}
