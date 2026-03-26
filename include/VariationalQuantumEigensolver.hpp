/*
 * Project: qubby
 * File name: VariationalQuantumEigensolver.hpp
 * Author: Cannelle Gourdet - lankley
 * File description: Declares the VQE class, which uses a parametric ansatz and the
 *                   parameter-shift gradient to minimise a Pauli Hamiltonian's energy.
 */

#ifndef VARIATIONALQUANTUMEIGENSOLVER_HPP_
    #define VARIATIONALQUANTUMEIGENSOLVER_HPP_

    #include "QuantumState.hpp"
    #include <vector>
    #include <functional>

/**
 * @enum PauliOp
 * @brief Single-qubit Pauli operator used to build Hamiltonian terms.
 *
 * - I: identity (no operation on this qubit).
 * - X, Y, Z: Pauli-X, Y, Z operators.
 */
enum class PauliOp { I, X, Y, Z };

/**
 * @struct PauliTerm
 * @brief One term of a Pauli Hamiltonian: a scalar coefficient times a tensor product of Pauli operators.
 *
 * Example — the Ising term -0.5 · Z₀Z₁:
 * @code
 *   PauliTerm { -0.5, { PauliOp::Z, PauliOp::Z } }
 * @endcode
 *
 * @var coeff Real coefficient of the term.
 * @var ops   Per-qubit Pauli operator (length must equal n_qubits).
 */
struct PauliTerm {
    double coeff;
    std::vector<PauliOp> ops;
};

/**
 * @class VQE
 * @brief Implements the Variational Quantum Eigensolver (VQE).
 *
 * VQE is a hybrid classical/quantum algorithm that finds the ground-state energy
 * of a Hamiltonian expressed as a sum of Pauli terms:
 * @code
 *   H = Σᵢ cᵢ · P₀ ⊗ P₁ ⊗ ... ⊗ Pₙ₋₁
 * @endcode
 *
 * A parametric ansatz |ψ(θ)⟩ is prepared on the quantum simulator, the energy
 * ⟨ψ(θ)|H|ψ(θ)⟩ is estimated by summing Pauli expectations, and θ is updated
 * classically via gradient descent using the **parameter-shift rule**:
 * @code
 *   ∂E/∂θᵢ = [E(θᵢ + π/2) - E(θᵢ - π/2)] / 2
 * @endcode
 *
 * Ansatz: Ry(θᵢ) on each qubit, then a linear chain of CNOT gates.
 * Optimiser: gradient descent, 200 iterations, learning rate 0.1.
 */
class VQE {
    private:
        /** @brief Number of qubits in the ansatz and Hamiltonian. */
        int n_qubits;

        /** @brief The Hamiltonian to minimise, as a list of Pauli terms. */
        std::vector<PauliTerm> hamiltonian;

        /** @brief Current best parameters θ (updated by run()). */
        std::vector<double> params;

        /**
         * @brief Prepares the ansatz state |ψ(θ)⟩ on @p state.
         *
         * Applies Ry(θᵢ) to each qubit i, then CNOT(i → i+1) for i in [0, n-2].
         *
         * @param state QuantumState to initialise and modify in-place.
         * @param theta Parameter vector (length n_qubits).
         */
        void prepareAnsatz(QuantumState &state, const std::vector<double> &theta) const;

        /**
         * @brief Returns ⟨ψ(θ)|H|ψ(θ)⟩ — the total energy expectation.
         *
         * Sums measurePauliTerm() over all terms in @p hamiltonian.
         *
         * @param theta Parameter vector.
         * @return Energy estimate as a double.
         */
        double expectationValue(const std::vector<double> &theta) const;

        /**
         * @brief Returns the expectation value of a single Pauli term.
         *
         * Prepares |ψ(θ)⟩, rotates into the measurement basis (H for X, Rx(π/2) for Y),
         * then evaluates ⟨P⟩ analytically from the amplitude vector:
         * ⟨P⟩ = Σₓ |ψ(x)|² · sign(x), where sign(x) = (−1)^(parity of active bits).
         *
         * @param term Pauli term to measure.
         * @param theta Parameter vector.
         * @return coeff · ⟨P⟩.
         */
        double measurePauliTerm(const PauliTerm &term, const std::vector<double> &theta) const;

    public:
        /**
         * @brief Constructs a VQE instance for the given Hamiltonian.
         *
         * @param n_qubits    Number of qubits.
         * @param hamiltonian List of Pauli terms defining H.
         */
        VQE(int n_qubits, std::vector<PauliTerm> hamiltonian);

        /**
         * @brief Runs the VQE optimisation loop and returns the estimated ground-state energy.
         *
         * Initialises θ to 0.1 for all parameters, then runs the **Adam** optimiser
         * (β₁=0.9, β₂=0.999, ε=1e-8, lr=0.01) for up to 500 iterations.
         * Gradients are computed via the parameter-shift rule. Early stopping triggers
         * when the gradient L2 norm falls below 1e-6. Updates `params` in-place.
         *
         * @return Minimum energy ⟨ψ(θ*)|H|ψ(θ*)⟩ found.
         */
        double run();

        /**
         * @brief Returns the optimal parameter vector found by run().
         *
         * @return Reference to the internal `params` vector after optimisation.
         */
        std::vector<double> getOptimalParams() const;
};

#endif /* VARIATIONALQUANTUMEIGENSOLVER_HPP_ */
