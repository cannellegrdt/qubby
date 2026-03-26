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
#include <random>
#include <chrono>
#include <algorithm>
#include <omp.h>
#include <iostream>
#include <iomanip>

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

    #pragma omp parallel for
    for (long long i=0; i<loop_size; i++) {
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

void QuantumState::rxGate(double theta, int targetQubit) {
    double val_cos = cos(theta / 2);
    std::complex<double> val_sin = std::complex<double>(0, -1) * sin(theta / 2);
    Matrix quantumGate = {val_cos, val_sin, val_sin, val_cos};
    applyGate(targetQubit, quantumGate);
}

void QuantumState::ryGate(double theta, int targetQubit) {
    double val_cos = cos(theta / 2);
    double val_sin = sin(theta / 2);
    Matrix quantumGate = {val_cos, -val_sin, val_sin, val_cos};
    applyGate(targetQubit, quantumGate);
}

void QuantumState::rzGate(double theta, int targetQubit) {
    std::complex<double> val_neg_exp = std::exp(std::complex<double>(0, -1) * theta / 2.0);
    std::complex<double> val_pos_exp = std::exp(std::complex<double>(0, 1) * theta / 2.0);
    Matrix quantumGate = {val_neg_exp, 0, 0, val_pos_exp};
    applyGate(targetQubit, quantumGate);
}

/**
 * @details
 * Reuses the same bitwise index reconstruction as applyGate(): for each half-index i,
 * i0 and i1 are the two basis states that differ only at @p targetQubit.
 * The control condition is checked on i0 (bit @p controlQubit of i0).
 * Because i0 and i1 only differ at @p targetQubit (and @p controlQubit ≠ @p targetQubit),
 * the control bit has the same value in both - checking i0 is sufficient.
 * When the control bit is set, the two amplitudes are swapped (equivalent to applying X).
 */
void QuantumState::cnotGate(int controlQubit, int targetQubit) {
    long long loop_size = (1LL << (this->num_qubits-1));

    for (long long i=0; i<loop_size; i++) {
        uint64_t mask = (1ULL << targetQubit) - 1;
        uint64_t i0 = ((i >> targetQubit) << (targetQubit + 1)) | (i & mask);
        uint64_t i1 = i0 | (1ULL << targetQubit);

        bool isControlBitSet = (i0 & (1ULL << controlQubit)) != 0;
        if (isControlBitSet) {
            std::complex<double> temp_amp_i0 = amplitudes[i0];
            amplitudes[i0] = amplitudes[i1];
            amplitudes[i1] = temp_amp_i0;
        }
    }
}

/// @details CNOT-decomposition: CNOT(q0→q1), CNOT(q1→q0), CNOT(q0→q1).
///          Each CNOT call runs the full loop, so total cost is 3 × O(2^(n-1)).
void QuantumState::swapGate(int q0, int q1) {
    cnotGate(q0, q1);
    cnotGate(q1, q0);
    cnotGate(q0, q1);
}

/**
 * @details
 * Same bitwise index reconstruction as cnotGate(): each half-index i yields i0/i1
 * (basis states differing only at @p targetQubit). The flip condition requires
 * both control bits to be 1 in i0:
 * @code
 *   shouldFlip = (i0 & (1 << c0)) && (i0 & (1 << c1))
 * @endcode
 * As with cnotGate(), checking i0 is sufficient because i0 and i1 only differ
 * at @p targetQubit (c0 ≠ targetQubit and c1 ≠ targetQubit).
 */
void QuantumState::toffoliGate(int c0, int c1, int targetQubit) {
    long long loop_size = (1LL << (this->num_qubits-1));

    for (long long i=0; i<loop_size; i++) {
        uint64_t mask = (1ULL << targetQubit) - 1;
        uint64_t i0 = ((i >> targetQubit) << (targetQubit + 1)) | (i & mask);
        uint64_t i1 = i0 | (1ULL << targetQubit);

        bool shouldFlip = (i0 & (1ULL << c0)) && (i0 & (1ULL << c1));
        if (shouldFlip) {
            std::complex<double> temp_amp_i0 = amplitudes[i0];
            amplitudes[i0] = amplitudes[i1];
            amplitudes[i1] = temp_amp_i0;
        }
    }
}

/**
 * @details
 * 1. Build a probability vector: p[i] = |ψ[i]|² = std::norm(amplitudes[i]).
 * 2. Draw one sample from std::discrete_distribution (non-uniform, weighted by p).
 * 3. Collapse the state: zero all amplitudes, then set amplitudes[result] = 1.0.
 *
 * A static std::mt19937 is seeded once at startup (random_device XOR'd with the
 * current time) so that each call advances the same generator rather than
 * re-seeding it, which would produce identical sequences on platforms where
 * std::random_device is deterministic.
 *
 * @note The state is mutated after this call - it is no longer a superposition.
 *       Call initialize() to reset, or save a copy before measuring.
 */
int QuantumState::measure() {
    std::vector<double> probabilities;
    for (const auto &amp : amplitudes) {
        probabilities.push_back(std::norm(amp));
    }
    static std::mt19937 gen = []() {
        std::random_device rd;
        auto seed = rd() ^ static_cast<uint32_t>(
            std::chrono::high_resolution_clock::now().time_since_epoch().count());
        return std::mt19937(seed);
    }();
    std::discrete_distribution<> dist(probabilities.begin(), probabilities.end());
    int index = dist(gen);

    std::fill(amplitudes.begin(), amplitudes.end(), std::complex<double>(0.0, 0.0));
    amplitudes[index] = 1.0;
    
    return index;
}

/**
 * @details
 * 1. Accumulate P(k=1) = Σ |ψ[i]|² over all i where bit k is 1.
 * 2. Draw from Bernoulli(P(k=1)) to get the outcome (0 or 1).
 * 3. Compute the renormalisation factor 1/√P(outcome).
 * 4. Zero all amplitudes inconsistent with the outcome; scale the rest.
 *
 * The remaining state is a valid normalised quantum state conditioned on the
 * measurement outcome (partial collapse).
 */
int QuantumState::measureQubit(int k) {
    double p1 = 0.0;
    int num_states = amplitudes.size();

    for (int i = 0; i < num_states; ++i) {
        if ((i >> k) & 1)
            p1 += std::norm(amplitudes[i]);
    }

    static std::mt19937 gen = []() {
        std::random_device rd;
        auto seed = rd() ^ static_cast<uint32_t>(
            std::chrono::high_resolution_clock::now().time_since_epoch().count());
        return std::mt19937(seed);
    }();
    std::bernoulli_distribution dist(p1);
    int index = dist(gen) ? 1 : 0;

    double p_result = (index == 1) ? p1 : (1.0 - p1);
    double normalization_factor = 1.0 / std::sqrt(p_result);

    for (int i=0; i<num_states; i++) {
        int bit_k = (i >> k) & 1;
        if (bit_k == index) {
            amplitudes[i] *= normalization_factor;
        } else {
            amplitudes[i] = 0.0;
        }
    }
    return index;
}

/// @details Single-line multiply: amplitudes[index] *= -1.0. O(1).
void QuantumState::phaseFlip(int index) {
    amplitudes[index] *= -1.0;
}

/**
 * @details
 * Applies a -1 phase to all basis states except |0⟩ by negating amplitudes[1..2^n-1].
 * This is the second step of Grover's diffusion operator (between two H layers).
 */
void QuantumState::phaseFlipAllExceptZero() {
    for (size_t i = 1; i < amplitudes.size(); i++)
        amplitudes[i] *= -1.0;
}

/**
 * @details
 * Uses the same bitwise index reconstruction as applyGate(): iterates over 2^(n-1)
 * half-indices, reconstructs i0/i1 relative to @p target, then applies the phase
 * e^(iθ) to i1 only when the @p control bit in i0 is 1.
 *
 * Only i1 (the |1⟩ component of @p target) is modified because Rz leaves |0⟩ unchanged
 * up to a global phase. The loop is parallelised with OpenMP.
 */
void QuantumState::applyControlledRz(int control, int target, double theta) {
    long long loop_size = (1LL << (this->num_qubits-1));

    #pragma omp parallel for
    for (long long i=0; i<loop_size; i++) {
        uint64_t mask = (1ULL << target) - 1;
        uint64_t i0 = ((i >> target) << (target + 1)) | (i & mask);
        uint64_t i1 = i0 | (1ULL << target);

        if ((i0 >> control) & 1)
            amplitudes[i1] *= std::exp(std::complex<double>(0, theta));
    }
}

std::complex<double> QuantumState::getAmplitude(int i) const {
    return amplitudes[i];
}

/**
 * @details
 * Iterates over all 2^n amplitudes. For each i where std::norm(amplitudes[i]) > 1e-6:
 * - Converts i to a big-endian binary string (bit num_qubits-1 first).
 * - Formats the complex amplitude as `real+imagi` (3 decimal places).
 * - Prints probability as a percentage.
 * Basis states with negligible probability (numerical noise) are silently skipped.
 */
void QuantumState::printState() {
    int loop_size = amplitudes.size();
    for (int i=0; i<loop_size; i++) {
        double probability = std::norm(amplitudes[i]);

        if (probability > 1e-6) {
            std::string binary = "";
            for (int b = num_qubits-1; b>=0; b--) {
                binary += ((i >> b) & 1) ? '1' : '0';
            }

            std::cout << "|" << std::setw(num_qubits) << i << "⟩"
                    << "  " << binary
                    << "  " << std::fixed << std::setprecision(3) << amplitudes[i].real() << (amplitudes[i].imag() >= 0 ? "+" : "") << amplitudes[i].imag() << "i"
                    << "  (" << (probability * 100.0) << "%)\n";
        }
    }
}

void QuantumState::setAmplitude(int i, std::complex<double> val) {
    amplitudes[i] = val;
}
