/*
 * Project: qubby
 * File name: QuantumState.cpp
 * Author: Cannelle Gourdet - lankley
 * File description: Implements the QuantumState methods.
 */

#include "QuantumState.hpp"
#include <stdexcept>
#include <string>

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
