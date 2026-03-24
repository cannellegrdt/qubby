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
         * @brief Returns the complex amplitude at the given basis-state index.
         *
         * @param index Basis-state index in [0, 2^n).
         * @return The complex amplitude ψ[index].
         */
        std::complex<double> get(int index);
};

#endif /* QUANTUMSTATE_HPP_ */
