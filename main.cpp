#include <iostream>
#include <cmath>
#include <chrono>
#include <vector>
// Карпов: carpson@mail.ru

using namespace std;
using namespace std::chrono;

double pi_approx_simple(long N) {
    double answer = 0.5 / (double) N;
    double n = 1. / (double) N, nn = n * n;
    double isq = 0;
    for (long i = 1; i < N; ++i) {
        isq += i + i - 1;
        answer += std::sqrt(1 - isq * nn); // Optimize, motherfucker!
    }
    return answer * 4 * n;
}

int main() {
    std::cout.precision(17);
    std::vector<int> times;
    int N = 10;
    for (int i = 0; i < N; ++i) {
        auto time = std::chrono::steady_clock::now();
        std::cout << pi_approx_simple(1000000000) << std::endl;
        times.push_back((duration_cast<milliseconds>(steady_clock::now() - time)).count());
    }
    double mean = 0, error = 0;
    for (int t : times)
        mean += t;
    mean /= N;
    for (int t : times)
        error += (t - mean) * (t - mean);
    error = std::sqrt(error / N);
    std::cout << "Profiling results:\nRun time " << mean << " ± " << error << " ms, " << N << " samples" << std::endl;
    return 0;
}
