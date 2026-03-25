/*
 * Project: qubby
 * File name: DeutschJozsa.cpp
 * Author: Cannelle Gourdet - lankley
 * File description: Implements the Deutsch-Jozsa algorithm: one oracle call suffices
 *                   to determine whether f:{0,1}^n → {0,1} is constant or balanced.
 */

#include "DeutschJozsa.hpp"

DeutschJozsa::DeutschJozsa(int n, std::function<void(QuantumState&)> oracle)
    : n(n), oracle(oracle) {}

/**
 * @details
 * Circuit execution (n+1 qubits, qubit n is the ancilla):
 * 1. Initialize: |0⟩^⊗(n+1).
 * 2. Flip ancilla: X on qubit n → |0⟩^⊗n |1⟩.
 * 3. H on all n+1 qubits → equal superposition with ancilla in |−⟩.
 * 4. Apply oracle(state): encodes f as (−1)^f(x)|x⟩|−⟩ (phase kickback).
 * 5. H on input qubits 0..n-1 (not ancilla).
 * 6. Measure each input qubit:
 *    - All 0 → CONSTANT (quantum interference cancelled all non-zero terms).
 *    - Any 1 → BALANCED.
 *
 * @note measureQubit() is used instead of measure() so the ancilla qubit is
 *       left untouched and only the n input qubits are sampled.
 */
FunctionType DeutschJozsa::run() {
    state.initialize(n+1);
    state.xGate(n);

    for (int i=0; i<=n; i++)
        state.hGate(i);

    oracle(state);

    for (int i=0; i<n; i++)
        state.hGate(i);

    for (int i=0; i<n; i++) {
        if (state.measureQubit(i) == 1) return FunctionType::BALANCED;
    }
    return FunctionType::CONSTANT;
}
