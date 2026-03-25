/*
 * Project: qubby
 * File name: benchmarks/benchmark.cpp
 * Author: Cannelle Gourdet - lankley
 * File description: Wall-clock benchmark comparing applyGate (OpenMP) vs
 *                   applyGateSYCL (Intel Iris Xe) on qubit counts 20–25.
 *                   Compile with -DSYCL_ENABLED to include the SYCL column.
 */

#include "QuantumState.hpp"
#include "Matrix.hpp"

#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

using Clock = std::chrono::steady_clock;
using Ms = std::chrono::duration<double, std::milli>;

static const Matrix HADAMARD = {
    1.0 / std::sqrt(2.0),  1.0 / std::sqrt(2.0),
    1.0 / std::sqrt(2.0), -1.0 / std::sqrt(2.0)
};

/**
 * @brief Measures the minimum wall-clock time (ms) for calling fn() `reps` times.
 *
 * The state is re-initialised before each repetition to keep the benchmark
 * deterministic and avoid accumulating state drift.
 */
template <typename Fn>
static double measure_min_ms(int n_qubits, int reps, Fn fn) {
    double best = std::numeric_limits<double>::max();
    for (int r = 0; r < reps; ++r) {
        QuantumState state;
        state.initialize(n_qubits);

        auto t0 = Clock::now();
        fn(state);
        auto t1 = Clock::now();

        double elapsed = Ms(t1 - t0).count();
        if (elapsed < best)
            best = elapsed;
    }
    return best;
}

int main() {
    constexpr int QUBIT_MIN = 20;
    constexpr int QUBIT_MAX = 25;
    constexpr int REPS = 5;

    std::cout << "\n";
    std::cout << "  qubby - applyGate benchmark  (best of " << REPS << " runs)\n";
    std::cout << "  Gate: Hadamard on qubit 0\n\n";

#ifdef SYCL_ENABLED
    std::cout << "  " << std::setw(7)  << "qubits"
              << "  " << std::setw(12) << "states"
              << "  " << std::setw(14) << "OpenMP (ms)"
              << "  " << std::setw(14) << "SYCL (ms)"
              << "  " << std::setw(10) << "speedup"
              << "\n";
    std::cout << "  " << std::string(7,  '-')
              << "  " << std::string(12, '-')
              << "  " << std::string(14, '-')
              << "  " << std::string(14, '-')
              << "  " << std::string(10, '-')
              << "\n";
#else
    std::cout << "  " << std::setw(7)  << "qubits"
              << "  " << std::setw(12) << "states"
              << "  " << std::setw(14) << "OpenMP (ms)"
              << "\n";
    std::cout << "  " << std::string(7,  '-')
              << "  " << std::string(12, '-')
              << "  " << std::string(14, '-')
              << "\n";
#endif

    for (int n = QUBIT_MIN; n <= QUBIT_MAX; ++n) {
        long long num_states = 1LL << n;

        double omp_ms = measure_min_ms(n, REPS, [&](QuantumState& s) {
            s.applyGate(0, HADAMARD);
        });

#ifdef SYCL_ENABLED
        double sycl_ms = measure_min_ms(n, REPS, [&](QuantumState& s) {
            s.applyGateSYCL(0, HADAMARD);
        });
        double speedup = omp_ms / sycl_ms;

        std::cout << "  " << std::setw(7)  << n
                  << "  " << std::setw(12) << num_states
                  << "  " << std::setw(14) << std::fixed << std::setprecision(3) << omp_ms
                  << "  " << std::setw(14) << std::fixed << std::setprecision(3) << sycl_ms
                  << "  " << std::setw(10) << std::fixed << std::setprecision(2) << speedup << "x"
                  << "\n";
#else
        std::cout << "  " << std::setw(7)  << n
                  << "  " << std::setw(12) << num_states
                  << "  " << std::setw(14) << std::fixed << std::setprecision(3) << omp_ms
                  << "\n";
#endif
    }

    std::cout << "\n";

#ifndef SYCL_ENABLED
    std::cout << "  [SYCL column disabled - rebuild with 'make benchmark_sycl' to compare]\n\n";
#endif

    return 0;
}
