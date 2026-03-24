/*
 * Project: qubby
 * File name: QuantumState.cpp
 * Author: Cannelle Gourdet - lankley
 * File description: Implements the QuantumState methods.
 */

#include "QuantumState.hpp"
#include <stdexcept>
#include <string>
#include <stdint.h>

/**
 * @details
 * Memory layout: each element is a `std::complex<double>` (16 bytes), so the
 * total footprint is 2^n × 16 bytes. The hard cap is MB_LIMIT (512 MB),
 * which corresponds to n = 25 qubits.
 *
 * After a successful call the state satisfies:
 *   - amplitudes.size() == 2^n
 *   - amplitudes[0] == 1.0 + 0i   (the |00...0⟩ basis state has unit probability)
 *   - amplitudes[k] == 0.0 + 0i   for k > 0
 */
void QuantumState::initialize(int n) {
    this->num_qubits = n;
    size_t size = 1ULL << n;

    size_t mb = size * sizeof(std::complex<double>);
    if (mb > MB_LIMIT)
        throw std::runtime_error("Requested memory (" + std::to_string(mb) + "MB) exceeds limit (512 MB)");

    amplitudes.resize(size);
    amplitudes[0] = 1.0;
}

/**
 * @details
 * The loop runs over 2^(n-1) "half-indices" (all bit patterns except the targetQubit bit).
 * For each half-index i, the full basis-state indices are reconstructed via bitwise
 * insertion of a 0 or 1 at position targetQubit:
 *   - mask = (1 << targetQubit) - 1     → preserves lower bits
 *   - i0   = insert 0 at bit targetQubit
 *   - i1   = i0 | (1 << targetQubit)   → insert 1 at bit targetQubit
 *
 * Both amplitudes are read before either is written to avoid overwriting a value
 * that is still needed.
 */
void QuantumState::applyGate(int targetQubit, Matrix quantumGate) {
    long long loop_size = (1LL << (this->num_qubits-1));

    for (int i=0; i<loop_size; i++) {
        uint64_t mask = (1ULL << targetQubit) - 1;
        uint64_t i0 = ((i >> targetQubit) << (targetQubit + 1)) | (i & mask);
        uint64_t i1 = i0 | (1ULL << targetQubit);

        std::complex<double> new_i0 = quantumGate.m00 * amplitudes[i0] + quantumGate.m01 * amplitudes[i1];
        std::complex<double> new_i1 = quantumGate.m10 * amplitudes[i0] + quantumGate.m11 * amplitudes[i1];

        amplitudes[i0] = new_i0;
        amplitudes[i1] = new_i1;
    }
}

/// @details Gate matrix: X = { {0,1}, {1,0} }.
void QuantumState::xGate(int targetQubit) {
    Matrix quantumGate = {0, 1, 1, 0};
    applyGate(targetQubit, quantumGate);
}

/// @details Gate matrix: H = 1/√2 · { {1,1}, {1,-1} }. Uses `sqrt(2)` from `<cmath>`.
void QuantumState::hGate(int targetQubit) {
    double val = 1 / sqrt(2);
    Matrix quantumGate = {val, val, val, -val};
    applyGate(targetQubit, quantumGate);
}

/// @details Gate matrix: Z = { {1,0}, {0,-1} }.
void QuantumState::zGate(int targetQubit) {
    Matrix quantumGate = {1, 0, 0, -1};
    applyGate(targetQubit, quantumGate);
}
