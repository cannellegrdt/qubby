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
         *              | 1  -1 |
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

        /**
         * @brief Applies the Rx(θ) rotation gate around the X-axis of the Bloch sphere.
         *
         * @code
         *   Rx(θ) = | cos(θ/2)    -i·sin(θ/2) |
         *           | -i·sin(θ/2)  cos(θ/2)   |
         * @endcode
         *
         * @param theta    Rotation angle in radians.
         * @param targetQubit Index of the qubit to rotate (0-based).
         */
        void rxGate(double theta, int targetQubit);

        /**
         * @brief Applies the Ry(θ) rotation gate around the Y-axis of the Bloch sphere.
         *
         * @code
         *   Ry(θ) = | cos(θ/2)  -sin(θ/2) |
         *           | sin(θ/2)   cos(θ/2) |
         * @endcode
         *
         * @param theta    Rotation angle in radians.
         * @param targetQubit Index of the qubit to rotate (0-based).
         */
        void ryGate(double theta, int targetQubit);

        /**
         * @brief Applies the Rz(θ) rotation gate around the Z-axis of the Bloch sphere.
         *
         * @code
         *   Rz(θ) = | e^(-iθ/2)  0          |
         *           | 0           e^(+iθ/2) |
         * @endcode
         *
         * @param theta    Rotation angle in radians.
         * @param targetQubit Index of the qubit to rotate (0-based).
         */
        void rzGate(double theta, int targetQubit);

        /**
         * @brief Applies the Controlled-NOT (CNOT) gate.
         *
         * Flips @p targetQubit if and only if @p controlQubit is |1⟩.
         * This is the fundamental two-qubit entangling gate.
         *
         * Truth table (control | target → target'):
         * @code
         *   0 | 0 → 0
         *   0 | 1 → 1
         *   1 | 0 → 1   (flipped)
         *   1 | 1 → 0   (flipped)
         * @endcode
         *
         * @param controlQubit Index of the control qubit (0-based).
         * @param targetQubit  Index of the qubit to flip conditionally (0-based).
         */
        void cnotGate(int controlQubit, int targetQubit);

        /**
         * @brief Swaps the quantum states of two qubits.
         *
         * Exchanges the amplitudes of @p q0 and @p q1 across all basis states.
         * Implemented as three successive CNOT gates (the standard CNOT-decomposition
         * of SWAP), so no additional bitwise kernel is needed:
         * @code
         *   SWAP(q0, q1) = CNOT(q0→q1) · CNOT(q1→q0) · CNOT(q0→q1)
         * @endcode
         *
         * @param q0 Index of the first qubit (0-based).
         * @param q1 Index of the second qubit (0-based).
         */
        void swapGate(int q0, int q1);

        /**
         * @brief Applies the Toffoli (CCX / Controlled-Controlled-NOT) gate.
         *
         * Flips @p targetQubit if and only if both @p c0 and @p c1 are |1⟩.
         * This is the universal reversible classical gate: any classical boolean
         * circuit can be built from Toffoli gates alone.
         *
         * Truth table (c0, c1 | target → target'):
         * @code
         *   0, 0 | * → *   (unchanged)
         *   0, 1 | * → *   (unchanged)
         *   1, 0 | * → *   (unchanged)
         *   1, 1 | 0 → 1   (flipped)
         *   1, 1 | 1 → 0   (flipped)
         * @endcode
         *
         * @param c0          Index of the first control qubit (0-based).
         * @param c1          Index of the second control qubit (0-based).
         * @param targetQubit Index of the qubit to flip conditionally (0-based).
         */
        void toffoliGate(int c0, int c1, int targetQubit);

        /**
         * @brief Measures the quantum state and collapses it to a single basis state.
         *
         * Samples a basis-state index according to the Born rule (probability = |ψ[i]|²),
         * then projects the state onto that outcome: all amplitudes are set to 0 except
         * the measured one, which is set to 1.0 (renormalisation after collapse).
         *
         * @return The measured basis-state index in [0, 2^n). The binary representation
         *         of this index encodes the classical bit value of each qubit
         *         (bit k == 1 means qubit k was measured as |1⟩).
         */
        int measure();

        /**
         * @brief Measures a single qubit and collapses the state accordingly.
         *
         * Computes P(qubit k = |1⟩) = Σ |ψ[i]|² for all i where bit k is 1,
         * draws a Bernoulli sample, then renormalises the remaining amplitudes
         * (zeroing all basis states inconsistent with the outcome).
         *
         * Unlike measure(), the other qubits remain in superposition.
         *
         * @param k Index of the qubit to measure (0-based).
         * @return 0 or 1 — the classical outcome of the measurement.
         */
        int measureQubit(int k);

        /**
         * @brief Flips the sign of a single basis-state amplitude (phase kickback).
         *
         * Equivalent to applying a Z-like phase of -1 to basis state @p index:
         * ψ[index] *= -1. Used as the oracle step in Grover's algorithm.
         *
         * @param index Basis-state index whose amplitude is negated.
         */
        void phaseFlip(int index);

        /**
         * @brief Flips the sign of all amplitudes except ψ[0] (inversion about the mean).
         *
         * Implements the diffusion operator's phase-flip step:
         * sets all amplitudes to 0 except ψ[0], which retains its original value.
         * Applied in Grover's diffusion between two layers of Hadamard gates.
         *
         * @note This is not a standard unitary gate — it is a partial projection
         *       used internally by Grover::applyDiffusion().
         */
        void phaseFlipAllExceptZero();

        /**
         * @brief Applies a controlled-Rz(θ) gate: rotates @p target by e^(iθ) if @p control is |1⟩.
         *
         * For each basis state where the @p control qubit is |1⟩, multiplies the
         * |1⟩ amplitude of @p target by e^(iθ). Used by QFT to apply the
         * conditional phase rotations Rk = Rz(2π/2^k).
         *
         * @param control Index of the control qubit (0-based).
         * @param target  Index of the target qubit to phase-rotate (0-based).
         * @param theta   Rotation angle in radians.
         */
        void applyControlledRz(int control, int target, double theta);

        #ifdef SYCL_ENABLED
            /**
             * @brief GPU-accelerated variant of applyGate() using SYCL.
             *
             * Offloads the gate application loop to the device selected by
             * `sycl::default_selector_v` (Intel Iris Xe via oneAPI on this machine).
             * The amplitude vector and the 2×2 gate are transferred to the device
             * via SYCL buffers, processed in parallel, then implicitly copied back
             * when the buffers go out of scope.
             *
             * The kernel applies the same bitwise index trick as applyGate():
             * each work-item handles one (i0, i1) pair independently,
             * which makes the loop embarrassingly parallel.
             *
             * @param targetQubit Index of the qubit to transform (0-based).
             * @param quantumGate The 2×2 unitary matrix to apply.
             *
             * @note Only available when compiled with `-DSYCL_ENABLED` and the
             *       Intel oneAPI Base Toolkit.
             */
            void applyGateSYCL(int targetQubit, Matrix quantumGate);
        #endif

        /**
         * @brief Returns the complex amplitude at basis-state index @p i.
         *
         * Intended for testing and inspection only.
         *
         * @param i Basis-state index in [0, 2^n).
         * @return The complex amplitude ψ[i].
         */
        std::complex<double> getAmplitude(int i) const;

        /**
         * @brief Prints the non-negligible amplitudes of the quantum state to stdout.
         *
         * For each basis state i with |ψ[i]|² > 1e-6, prints one line:
         * @code
         *   |<index>⟩  <binary>  <real>+<imag>i  (<probability>%)
         * @endcode
         * where `<binary>` is the big-endian bit string of i (MSB = qubit n-1).
         * Basis states with negligible probability are silently skipped.
         */
        void printState();

        void setAmplitude(int i, std::complex<double> val);
};

#endif /* QUANTUMSTATE_HPP_ */
