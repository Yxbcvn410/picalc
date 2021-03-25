#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include <vector>
#include "mpi/mpi.h"
// Карпов: carpson@mail.ru

using namespace std;
using namespace std::chrono;

double pi_approx_simple(long N, int k, int rk) {
    double answer = 0.5 / (double) N;
    double n = 1. / (double) N, nn = n * n;
    if (rk == 0)
        rk = k;
    for (long i = rk; i < N; i += k) {
        answer += std::sqrt(1 - (double) (i * i) * nn);
    }
    return answer * 4 * n;
}

int main() {
    int k, rk, N;
    MPI_Comm_size(MPI_COMM_WORLD, &k);
    MPI_Comm_rank(MPI_COMM_WORLD, &rk);

    if (rk == 0) {
        ifstream fs("N.dat", ios::in);
        fs >> N;
        fs.close();
        for (int i = 1; i < k; ++i)
            MPI_Send(&N, sizeof(typeof(N)), MPI_INT, i, 123, MPI_COMM_WORLD);
    } else {
        MPI_Recv(&N, sizeof(typeof(N)), MPI_INT, 0, 123, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    float ps = pi_approx_simple(N, k, rk);

    if (rk == 0){
        float msg;
        for (int i = 0; i < N; ++i) {
            MPI_Recv(&msg, sizeof(typeof(msg)), MPI_FLOAT, i, 234, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            ps += msg;
        }
        cout << msg;
    } else {
        MPI_Send(&ps, sizeof(typeof(ps)), MPI_FLOAT, 0, 234, MPI_COMM_WORLD);
    }
}
