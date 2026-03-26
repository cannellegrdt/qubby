/*
 * Project: qubby
 * File name: QuantumFourierTransform.cpp
 * Author: Cannelle Gourdet - lankley
 * File description: Implements the Quantum Fourier Transform via Hadamard + controlled-Rz
 *                   gates and bit-reversal swaps.
 */

#include "QuantumFourierTransform.hpp"
#include "QuantumState.hpp"

QFT::QFT(int n_qubits) {
    state.initialize(n_qubits);
    this->n_qubits = n_qubits;
}

/// @details Returns a mutable reference so callers can inspect or chain operations
///          on the transformed state without copying the amplitude vector.
QuantumState& QFT::getState() {
    return state;
}

/**
 * @details
 * The circuit follows the standard recursive QFT decomposition:
 *
 * For each qubit j from 0 to n-1:
 *   1. H(j): puts qubit j into superposition.
 *   2. For each k from j+1 to n-1:
 *      CRz(2π / 2^(k-j+1)) controlled by k, targeting j.
 *      This applies the conditional phase rotation Rk that encodes
 *      the fractional binary expansion of the input's Fourier coefficient.
 *
 * After the nested loops the qubit order is reversed relative to the standard
 * QFT definition, so n/2 SWAP gates restore the correct bit order:
 *   swap(0, n-1), swap(1, n-2), ...
 *
 * Total gate count: n·H + n(n-1)/2·CRz + floor(n/2)·SWAP = O(n²).
 */
void QFT::run() {
    applyQFT(this->state, this->n_qubits);
}

void QFT::applyQFT(QuantumState& s, int n_qubits) {
    for (int j=0; j<n_qubits; j++) {
        s.hGate(j);
        for (int k=j+1; k<n_qubits; k++) {
            double theta = (2.0 * M_PI) / (1LL << (k - j + 1));
            s.applyControlledRz(k, j, theta);
        }
    }
    int loop_size = n_qubits/2;
    for (int i=0; i<loop_size; i++) {
        s.swapGate(i, n_qubits-1-i);
    }
}
