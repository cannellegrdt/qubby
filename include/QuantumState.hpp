/*
 * Project: qubby
 * File name: QuantumState.hpp
 * Author: Cannelle Gourdet - lankley
 * File description: Declares the QuantumState class, which represents a multi-qubit quantum state
 *                   as a vector of complex amplitudes, with a 512 MB memory cap.
 */

#ifndef QUANTUMSTATE_HPP_
    #define QUANTUMSTATE_HPP_

    #include <vector>
    #include <complex>

    #include "Matrix.hpp"

/**
 * @brief Maximum allowed memory for the state vector (512 MB).
 *
 * Prevents accidental allocation of excessively large state vectors.
 * At 20 qubits the vector requires 16 MB; at 25 qubits it requires 512 MB.
 */
constexpr size_t MB_LIMIT = 512ULL * 1024 * 1024;

/**
 * @class QuantumState
 * @brief Represents the quantum state of an n-qubit system as a complex amplitude vector.
 *
 * The state is stored as a contiguous array of 2^n complex amplitudes.
 * The initial state is always |00...0⟩ (only amplitude[0] == 1.0, all others == 0.0).
 *
 * Gates are applied in-place via bitwise index manipulation to avoid
 * constructing full 2^n × 2^n matrices.
 *
 * @note Supports up to 25 qubits (512 MB). Attempting more will throw.
 */
class QuantumState {
    private:
        /** @brief Complex probability amplitudes of the quantum state (size = 2^n). */
        std::vector<std::complex<double>> amplitudes;

        /** @brief Number of qubits represented by this state. */
        int num_qubits;

    public:
        /**
         * @brief Initialises the quantum state to |00...0⟩ for @p n qubits.
         *
         * Allocates a vector of 2^n complex amplitudes, sets all to 0, then sets
         * amplitude[0] = 1.0 (the |00...0⟩ basis state). Throws if the required
         * memory would exceed MB_LIMIT.
         *
         * @param n Number of qubits (must satisfy 2^n * 16 bytes ≤ 512 MB).
         * @throws std::runtime_error if the required allocation exceeds MB_LIMIT.
         */
        void initialize(int n);

        /**
         * @brief Applies an arbitrary single-qubit gate to the target qubit.
         *
         * Uses the bitwise index trick to iterate over all 2^(n-1) pairs of basis
         * states that differ only in the @p targetQubit bit, and applies @p quantumGate
         * in-place on each pair. This avoids building the full 2^n × 2^n matrix.
         *
         * For a given loop index i, the two paired basis states are:
         *   - i0: index where bit @p targetQubit is 0
         *   - i1: index where bit @p targetQubit is 1 (i1 = i0 | (1 << targetQubit))
         *
         * The update is:
         * @code
         *   ψ[i0]' = m00·ψ[i0] + m01·ψ[i1]
         *   ψ[i1]' = m10·ψ[i0] + m11·ψ[i1]
         * @endcode
         *
         * @param targetQubit Index of the qubit to transform (0-based, 0 = least significant).
         * @param quantumGate The 2×2 unitary matrix to apply.
         */
        void applyGate(int targetQubit, Matrix quantumGate);

        /**
         * @brief Applies the Pauli-X (NOT) gate to the target qubit.
         *
         * Flips the |0⟩ and |1⟩ amplitudes of @p targetQubit:
         * @code
         *   X = | 0  1 |
         *       | 1  0 |
         * @endcode
         *
         * @param targetQubit Index of the qubit to flip (0-based).
         */
        void xGate(int targetQubit);

        /**
         * @brief Applies the Hadamard gate to the target qubit.
         *
         * Creates an equal superposition of |0⟩ and |1⟩:
         * @code
         *   H = 1/√2 · | 1   1 |
         *               | 1  -1 |
         * @endcode
         *
         * Applying H twice returns the qubit to its original state (H² = I).
         *
         * @param targetQubit Index of the qubit to put in superposition (0-based).
         */
        void hGate(int targetQubit);

        /**
         * @brief Applies the Pauli-Z (phase flip) gate to the target qubit.
         *
         * Leaves |0⟩ unchanged and flips the sign of |1⟩:
         * @code
         *   Z = | 1   0 |
         *       | 0  -1 |
         * @endcode
         *
         * @param targetQubit Index of the qubit to phase-flip (0-based).
         */
        void zGate(int targetQubit);
};

#endif /* QUANTUMSTATE_HPP_ */
