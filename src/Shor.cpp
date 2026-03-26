/*
 * Project: qubby
 * File name: Shor.cpp
 * Author: Cannelle Gourdet - lankley
 * File description: Implements Shor's quantum factorisation algorithm: order finding via QFT,
 *                   modular exponentiation oracle, and continued-fraction period extraction.
 */

#include "Shor.hpp"
#include "QuantumFourierTransform.hpp"
#include <stdint.h>
#include <vector>
#include <random>

/**
 * @brief Computes gcd(a, b) via Euclidean algorithm. O(log min(a,b)).
 */
long long gcd(long long a, long long b) {
    while (b != 0) {
        long long temp = b;
        b = a % b;
        a = temp;
    }
    return a;
}

/**
 * @brief Computes base^exp mod m via fast binary exponentiation. O(log exp).
 */
long long modpow(long long base, long long exp, long long mod) {
    long long result = 1;
    base = base % mod;

    while (exp > 0) {
        if (exp % 2 == 1) result = result * base % mod;
        base = base * base % mod;
        exp /= 2;
    }
    return result;
}

/**
 * @brief Extracts the period r from a QFT measurement m via continued fractions.
 *
 * Approximates m/Q as a continued fraction p/q and returns the first denominator q
 * satisfying |m/Q − p/q| < 1/(2Q) and q < N. This is the standard post-processing
 * step in Shor's order-finding subroutine.
 *
 * @param m Measured value from the QFT register (0 ≤ m < Q).
 * @param Q Size of the period-finding register (Q = 2^(2n)).
 * @param N Integer being factorised (used as upper bound on q).
 * @return Candidate period r, or 0 if no convergent meets the criteria.
 */
long long continuedFraction(long long m, long long Q, int N) {
    if (m == 0) return 0;

    double x = static_cast<double>(m) / Q;
    double target = 1.0 / (2.0 * Q);

    long long p_prev2 = 0, p_prev1 = 1;
    long long q_prev2 = 1, q_prev1 = 0;

    double remainder = x;

    while (true) {
        long long a = static_cast<long long>(std::floor(remainder));
        long long p_k = a * p_prev1 + p_prev2;
        long long q_k = a * q_prev1 + q_prev2;
        if (q_k >= N) break;

        double current_val = static_cast<double>(p_k) / q_k;
        if (std::abs(x-current_val) <= target) return q_k;

        p_prev2 = p_prev1;
        p_prev1 = p_k;
        q_prev2 = q_prev1;
        q_prev1 = q_k;

        double diff = remainder - a;
        if (std::abs(diff) < 1e-15) break;
        remainder = 1.0 / diff;
    }
    return 0;
}

/**
 * @details
 * Up to 10 attempts with a fresh random base a each time:
 * 1. If gcd(a, N) ≠ 1, a is already a factor — return immediately.
 * 2. Initialize 3n qubits; apply H^⊗2n to the period register.
 * 3. Apply modularOracle(a, N): encodes a^x mod N in the ancilla register.
 * 4. Measure ancilla qubits [2n, 3n-1]: collapses the period register to a
 *    uniform superposition over {x : a^x mod N == measured value}.
 * 5. Apply QFT to the period register via qft.applyQFT().
 * 6. Measure the full state; extract the lower 2n bits as m.
 * 7. Run continuedFraction(m, Q, N) to find r.
 * 8. If r is even and a^(r/2) ≢ -1 (mod N): try gcd(a^(r/2) ± 1, N).
 */
std::pair<long long, long long> Shor::run(int N) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<long long> distrib(2, N - 1);

    int n = 0;
    while ((1LL << n) < N) n++;
    long long Q = 1LL << (2*n);

    for (int attempt = 0; attempt < 10; attempt++) {
        state.initialize(3*n);
        for (int i=0; i<2*n; i++) state.hGate(i);

        long long a = distrib(gen);
        long long gcdAN = gcd(a, N);
        if (gcdAN != 1)
            return {gcdAN, N/gcdAN};
        modularOracle(a, N);

        for (int i=0; i<n; i++) state.measureQubit(2*n + i);

        QFT qft(2*n);
        qft.applyQFT(state, 2*n);

        long long measured = state.measure();
        long long m = measured & (Q - 1);
        long long r = continuedFraction(m, Q, N);
        if (r == 0) continue;
        if (r % 2 == 0 && modpow(a, r/2, N) != N-1) {
            long long p = gcd(modpow(a, r/2, N)-1, N);
            long long q = gcd(modpow(a, r/2, N)+1, N);
            if (p != 1 && p != N) return {p, q};
        }
    }
    return {-1, -1};
}

/**
 * @details
 * State layout: qubits [0, 2n-1] hold x (period register), qubits [2n, 3n-1] hold y (ancilla).
 * For each basis state |x⟩|y⟩ with non-negligible amplitude:
 *   - Extracts x = i & x_mask  (lower 2n bits).
 *   - Extracts y = i >> 2n     (upper n bits).
 *   - Computes fx = a^x mod N via modpow.
 *   - Maps to new_i = |x⟩|y ⊕ fx⟩ and accumulates the amplitude.
 * The result is written back via setAmplitude().
 */
void Shor::modularOracle(long long a, long long N) {
    int n = 0;
    while ((1LL << n) < N) n++;
    uint64_t total_states = 1ULL << (3*n);
    uint64_t x_mask = (1ULL << (2*n)) - 1;

    std::vector<std::complex<double>> new_amplitudes(total_states, 0.0);

    for (uint64_t i=0; i<total_states; i++) {
        std::complex<double> amp = state.getAmplitude(i);
        if (std::abs(amp) < 1e-15) continue;

        uint64_t x = i & x_mask;
        uint64_t y = i >> (2*n);

        uint64_t fx = modpow(a, x, N);
        uint64_t new_i = ((y ^ fx) << (2*n)) | x;
        new_amplitudes[new_i] += amp;
    }

    for (uint64_t i=0; i<total_states; i++)
        state.setAmplitude(i, new_amplitudes[i]);
}
