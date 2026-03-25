/*
 * Project: qubby
 * File name: Grover.cpp
 * Author: Cannelle Gourdet - lankley
 * File description: Implements Grover's amplitude-amplification algorithm,
 *                   searching a 2^n-element space in O(√(2^n)) oracle queries.
 */

#include "Grover.hpp"
#include "QuantumState.hpp"
#include <cmath>

Grover::Grover(int n_qubits, int target) {
    this->n_qubits = n_qubits;
    this->target = target;
}

/// @details Delegates to QuantumState::phaseFlip(target): ψ[target] *= -1. O(1).
void Grover::applyOracle() {
    state.phaseFlip(this->target);
}

/**
 * @details
 * The diffusion operator 2|s⟩⟨s| − I (reflection about the uniform superposition |s⟩)
 * is decomposed into three steps:
 * 1. H^⊗n: rotates from the computational basis to the Hadamard basis.
 * 2. phaseFlipAllExceptZero: applies a -1 phase to all states except |0⟩,
 *    equivalent to 2|0⟩⟨0| − I in the Hadamard basis.
 * 3. H^⊗n: rotates back to the computational basis.
 *
 * The net effect is an inversion of all amplitudes about their mean, which
 * doubles the target amplitude (roughly) each iteration.
 */
void Grover::applyDiffusion() {
    for (int i=0; i<this->n_qubits; i++)
        state.hGate(i);

    state.phaseFlipAllExceptZero();

    for (int i=0; i<this->n_qubits; i++)
        state.hGate(i);
}

/**
 * @details
 * 1. Initialize: |0⟩^⊗n.
 * 2. H^⊗n: equal superposition over all 2^n basis states.
 * 3. Iterate k = floor((π/4)·√(2^n)) times:
 *    - applyOracle()    → marks target with -1 phase
 *    - applyDiffusion() → amplifies target amplitude
 * 4. measure(): collapses state; returns target with probability ≥ 1 − 1/N.
 *
 * The number of iterations k is optimal: fewer or more iterations decrease
 * the success probability (the amplitude oscillates sinusoidally).
 */
int Grover::run() {
    state.initialize(n_qubits);

    for (int i=0; i<this->n_qubits; i++)
        state.hGate(i);

    int k = (int)(M_PI / 4.0 * std::sqrt((double)(1 << n_qubits)));

    for (int i=0; i<k; i++) {
        applyOracle();
        applyDiffusion();
    }
    return state.measure();
}
