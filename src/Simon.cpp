/*
 * Project: qubby
 * File name: Simon.cpp
 * Author: Cannelle Gourdet - lankley
 * File description: Implements Simon's algorithm: oracle application, repeated quantum sampling,
 *                   and GF(2) Gaussian elimination to recover the hidden period.
 */

#include "Simon.hpp"
#include <stdint.h>
#include <vector>

/**
 * @details
 * Gaussian elimination over GF(2) on the equation system {y · s = 0}:
 * 1. For each bit position k (0 to n-1), find a row with bit k set and record
 *    it as the pivot; XOR it into all other rows that have bit k set.
 * 2. After elimination, any bit position with no pivot is a free variable.
 *    The first free bit is set to 1 in s to produce a non-trivial solution.
 * 3. Back-substitution: for each pivot row, propagate the free-variable
 *    contribution into the remaining bits of s.
 *
 * Returns 0 (trivial solution s = 0) only if the system has full rank,
 * i.e., f is injective (no hidden period).
 */
int Simon::solveGF2(std::vector<int> &equations, int n) {
    int size_eq = equations.size();
    std::vector<int> pivot_row(n, -1);

    for (int k=0; k<n; k++) {
        int j_found = -1;

        for (int j = 0; j < size_eq; j++) {
            if (equations[j] & (1ULL << k)) {
                j_found = j;
                break;
            }
        }
        if (j_found != -1) {
            pivot_row[k] = j_found;

            for (int i=0; i<size_eq; i++) {
                if (i != j_found && (equations[i] & (1ULL << k)))
                    equations[i] ^= equations[j_found];
            }
        } 
    }

    uint64_t s = 0;

    for (int k = 0; k < n; k++) {
        if (pivot_row[k] == -1) {
            s |= (1ULL << k);
            break;
        }
    }

    for (int k=n-1; k>=0; k--) {
        if (pivot_row[k] == -1) continue;

        uint64_t row = equations[pivot_row[k]];
        for (int b = 0; b < n; b++) {
            if (b == k)
                continue;
            if ((row & (1ULL << b)) && (s & (1ULL << b)))
                s ^= (1ULL << k);
        }
    }
    return s;
}

/**
 * @details
 * Repeats the following until n-1 linearly independent equations are collected:
 * 1. initialize(2n): reset the 2n-qubit state.
 * 2. H^⊗n on the upper n qubits (input register).
 * 3. applyOracle(f, n): compute |x⟩|0⟩ → |x⟩|f(x)⟩.
 * 4. measureQubit() on the lower n qubits: collapses to a fixed f(x₀),
 *    leaving the input register in 1/√2 (|x₀⟩ + |x₀ ⊕ s⟩).
 * 5. H^⊗n on the upper n qubits.
 * 6. Measure each upper qubit to get y satisfying y · s = 0 (mod 2).
 * 7. Add y to the equation set only if it is linearly independent (Gaussian basis check).
 *
 * After n-1 independent equations, calls solveGF2() to recover s.
 */
int Simon::run(int n, std::function<int(int)> f) {
    if (n == 1) return (f(0) == f(1)) ? 1 : 0;

    std::vector<int> equations;
    std::vector<int> basis(n, 0);

    int max_iter = 100 * n;
    while ((int)equations.size() < n-1 && max_iter-- > 0) {
        state.initialize(2*n);

        for (int i=n; i<2*n; i++) state.hGate(i);
        applyOracle(f, n);
        for (int i=0; i<n; i++) state.measureQubit(i);
        for (int i=n; i<2*n; i++) state.hGate(i);

        int y = 0;
        for (int i=0; i<n; i++) {
            int k = state.measureQubit(n+i);
            y = y | (k << i);
        }

        int reduced = y;
        for (int b = n-1; b >= 0; b--) {
            if (!((reduced >> b) & 1)) continue;
            if (!basis[b]) {
                basis[b] = reduced;
                equations.push_back(y);
                break;
            }
            reduced ^= basis[b];
        }
    }
    return solveGF2(equations, n);
}

/**
 * @details
 * State layout: upper n bits = x (input), lower n bits = y (output).
 * For each basis state |x⟩|y⟩ with |amp| > 1e-12:
 *   - Extracts x = (i >> n) & mask.
 *   - Computes fx = f(x).
 *   - Maps to new_i = |x⟩|y ⊕ f(x)⟩ and accumulates amplitude.
 * The new amplitude vector is written back via setAmplitude().
 */
void Simon::applyOracle(std::function<int(int)> f, int n) {
    uint64_t loop_size = 1ULL << (2*n);
    std::vector<std::complex<double>> new_amplitudes(loop_size, 0.0);

    for (uint64_t i=0; i<loop_size; i++) {
        std::complex<double> amp = state.getAmplitude(i);
        if (std::abs(amp) < 1e-12) continue;

        int x = (i >> n) & ((1LL << n) - 1);
        int out = i & ((1LL << n) - 1);

        int fx = f(x);
        int new_out = out ^ fx;
        uint64_t new_i = ((uint64_t)x << n) | new_out;
        new_amplitudes[new_i] += amp;
    }

    for (uint64_t i=0; i<loop_size; i++)
        state.setAmplitude(i, new_amplitudes[i]);
}
