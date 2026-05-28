/*
 * Project: qubby
 * File name: DensityMatrix.hpp
 * Author: Cannelle Gourdet - lankley
 * File description: Declares the DensityMatrix class, which represents a noisy quantum state
 *                   as a 2^n × 2^n complex density matrix and supports Kraus-operator noise channels.
 */

#ifndef DENSITYMATRIX_HPP_
    #define DENSITYMATRIX_HPP_
    #include <vector>
    #include <complex>
    #include "Matrix.hpp"

/**
 * @brief Maximum number of qubits for density-matrix simulation.
 *
 * A density matrix for n qubits requires 4^n complex<double> values (16 bytes each):
 *   n=10 → 16 MB, n=11 → 64 MB, n=12 → 256 MB, n=13 → 1 GB.
 * Capped at 12 to stay within 512 MB.
 */
constexpr int DM_MAX_QUBITS = 12;

/**
 * @class DensityMatrix
 * @brief Represents a mixed quantum state via its density matrix ρ (size 2^n × 2^n).
 *
 * Unlike a pure-state vector, a density matrix can represent statistical mixtures
 * of quantum states produced by decoherence. Stored in row-major order:
 * element ρ_{ij} lives at rho[i * dim + j].
 *
 * Gates are applied as  ρ → U ρ U†  (left then right multiplication).
 * Noise channels are applied as  ρ → Σᵢ Kᵢ ρ Kᵢ†  (Kraus operator sum).
 *
 * Three noise channels are provided:
 *   - Bit-flip      with Kraus ops  {√(1-p)·I, √p·X}
 *   - Phase-flip    with Kraus ops  {√(1-p)·I, √p·Z}
 *   - Depolarizing  with Kraus ops  {√(1-p)·I, √(p/3)·X, √(p/3)·Y, √(p/3)·Z}
 */
class DensityMatrix {
    private:
        /** @brief Density matrix elements in row-major order (size = dim²). */
        std::vector<std::complex<double>> rho;

        /** @brief Number of qubits. */
        int num_qubits;

        /** @brief 2^num_qubits. */
        int dim;

        /**
         * @brief Applies U to qubit q on all row indices: ρ → (U_q ⊗ I) ρ.
         *
         * For each pair of row indices (i0, i1) differing only at bit q, and for
         * every column j, updates ρ[i0,j] and ρ[i1,j] using the 2×2 gate U.
         */
        void gateLeft(int q, Matrix U);

        /**
         * @brief Applies U† to qubit q on all column indices: ρ → ρ (U_q ⊗ I)†.
         *
         * For each pair of column indices (j0, j1) differing only at bit q, and for
         * every row i, updates ρ[i,j0] and ρ[i,j1] using the conjugate-transpose U†.
         */
        void gateRight(int q, Matrix U);

        /**
         * @brief Applies a list of Kraus operators to qubit q: ρ → Σᵢ Kᵢ ρ Kᵢ†.
         *
         * For each K in ops:
         *   1. Computes temp = K · ρ (left multiplication on rows).
         *   2. Accumulates new_ρ += temp · K† (right multiplication on columns).
         * Then sets ρ = new_ρ.
         */
        void kraus(int q, const std::vector<Matrix>& ops);

        /** @brief Left-multiplies by CNOT (control, target) on row indices. */
        void cnotLeft(int ctrl, int tgt);

        /** @brief Right-multiplies by CNOT† = CNOT on column indices. */
        void cnotRight(int ctrl, int tgt);

        /** @brief Left-multiplies by Toffoli on row indices. */
        void toffoliLeft(int c0, int c1, int tgt);

        /** @brief Right-multiplies by Toffoli† = Toffoli on column indices. */
        void toffoliRight(int c0, int c1, int tgt);

    public:
        /**
         * @brief Initialises the density matrix to the pure state |00...0⟩⟨00...0|.
         *
         * Allocates a dim×dim matrix (dim = 2^n), sets all elements to zero,
         * then sets ρ[0,0] = 1.
         *
         * @param n Number of qubits (must be ≤ DM_MAX_QUBITS).
         * @throws std::runtime_error if n > DM_MAX_QUBITS.
         */
        void initialize(int n);

        /**
         * @brief Applies a single-qubit gate U to qubit q: ρ → U_q ρ U_q†.
         */
        void applyGate(int q, Matrix U);

        void xGate(int q);
        void hGate(int q);
        void zGate(int q);
        void rxGate(double theta, int q);
        void ryGate(double theta, int q);
        void rzGate(double theta, int q);

        /**
         * @brief Applies CNOT gate: ρ → CNOT_{ctrl,tgt} ρ CNOT_{ctrl,tgt}†.
         */
        void cnotGate(int ctrl, int tgt);

        /** @brief Applies SWAP via three CNOTs. */
        void swapGate(int q0, int q1);

        /**
         * @brief Applies Toffoli gate: ρ → CCX_{c0,c1,tgt} ρ CCX_{c0,c1,tgt}†.
         */
        void toffoliGate(int c0, int c1, int tgt);

        /**
         * @brief Applies the bit-flip (X-error) channel to qubit q.
         *
         * Kraus operators: K0 = √(1-p)·I,  K1 = √p·X.
         * Models a physical bit-flip error occurring with probability p.
         *
         * @param q Qubit index.
         * @param p Error probability in [0, 1].
         */
        void applyBitFlip(int q, double p);

        /**
         * @brief Applies the phase-flip (Z-error) channel to qubit q.
         *
         * Kraus operators: K0 = √(1-p)·I,  K1 = √p·Z.
         * Models dephasing noise with probability p.
         *
         * @param q Qubit index.
         * @param p Error probability in [0, 1].
         */
        void applyPhaseFlip(int q, double p);

        /**
         * @brief Applies the depolarizing channel to qubit q.
         *
         * Kraus operators: K0 = √(1-p)·I,  K1 = √(p/3)·X,  K2 = √(p/3)·Y,  K3 = √(p/3)·Z.
         * Each Pauli error (X, Y, Z) occurs with probability p/3.
         * This is the most general single-qubit noise model and subsumes
         * bit-flip and phase-flip as special cases.
         *
         * @param q Qubit index.
         * @param p Total error probability in [0, 1].
         */
        void applyDepolarizing(int q, double p);

        /**
         * @brief Samples a measurement outcome from the diagonal of ρ and collapses the state.
         *
         * Probabilities are the diagonal elements ρ[i,i] (Born rule for mixed states).
         * After measurement, ρ = |result⟩⟨result| (pure-state projection).
         *
         * @return Measured basis-state index in [0, 2^n).
         */
        int measure();

        /**
         * @brief Prints basis-state probabilities from the diagonal of ρ.
         *
         * Only prints states with probability > 1e-6.
         */
        void printState();

        /**
         * @brief Returns the probability of basis state i: the i-th diagonal element ρ[i,i].
         *
         * Intended for testing and inspection only.
         *
         * @param i Basis-state index in [0, 2^n).
         * @return Real part of ρ[i,i] (always ≥ 0 for a valid density matrix).
         */
        double getDiagonal(int i) const;
};

#endif /* DENSITYMATRIX_HPP_ */
