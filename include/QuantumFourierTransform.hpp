/*
 * Project: qubby
 * File name: QuantumFourierTransform.hpp
 * Author: Cannelle Gourdet - lankley
 * File description: Declares the QFT class, which implements the Quantum Fourier Transform
 *                   as a sequence of Hadamard and controlled-Rz gates followed by SWAP reversals.
 */

#ifndef QUANTUMFOURIER_TRANSFORM_HPP_
    #define QUANTUMFOURIER_TRANSFORM_HPP_
    #include "QuantumState.hpp"

/**
 * @class QFT
 * @brief Implements the Quantum Fourier Transform on n qubits.
 *
 * The QFT is the quantum analogue of the Discrete Fourier Transform (DFT).
 * It maps basis state |j⟩ to:
 * @code
 *   QFT|j⟩ = (1/√N) Σ_{k=0}^{N-1} e^(2πijk/N) |k⟩ where N = 2^n
 * @endcode
 *
 * Circuit for qubit j:
 * @code
 *   H(j)  →  CRz(2π/2²) with j+1  →  CRz(2π/2³) with j+2  →  ...
 * @endcode
 * Followed by bit-reversal via SWAP gates.
 *
 * The QFT is a core building block of Shor's algorithm and quantum phase estimation.
 */
class QFT {
    private:
        /** @brief Internal quantum state of n qubits. */
        QuantumState state;

        /** @brief Number of qubits in the transform. */
        int n_qubits;

    public:
        /**
         * @brief Constructs a QFT instance and initialises the state to |00...0⟩.
         *
         * @param n_qubits Number of qubits to transform.
         * @throws std::runtime_error if the required memory exceeds MB_LIMIT.
         */
        QFT(int n_qubits);

        /**
         * @brief Applies the QFT circuit in-place on the internal state.
         *
         * For each qubit j (0 to n-1):
         *   1. Apply H(j).
         *   2. For each k > j: apply CRz(2π/2^(k-j+1)) controlled by k, targeting j.
         * Then reverse qubit order with n/2 SWAP gates.
         */
        void run();

        /**
         * @brief Returns a reference to the internal QuantumState after the transform.
         *
         * Allows inspection of amplitudes or chaining with other operations.
         *
         * @return Reference to the internal QuantumState.
         */
        QuantumState &getState();

        void applyQFT(QuantumState& s, int n_qubits);
};

#endif /* QUANTUMFOURIER_TRANSFORM_HPP_ */
