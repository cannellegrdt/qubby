/*
 * Project: qubby
 * File name: DeutschJozsa.hpp
 * Author: Cannelle Gourdet - lankley
 * File description: Declares the DeutschJozsa class, which determines in a single quantum
 *                   query whether a boolean function f:{0,1}^n → {0,1} is constant or balanced.
 */

#ifndef DEUTSCHJOZSA_HPP_
    #define DEUTSCHJOZSA_HPP_

    #include "QuantumState.hpp"
    #include <functional>

/**
 * @enum FunctionType
 * @brief Classification of a boolean oracle function.
 *
 * - CONSTANT: f(x) is the same for all inputs (all 0 or all 1).
 * - BALANCED: f(x) is 0 for exactly half the inputs and 1 for the other half.
 */
enum class FunctionType { CONSTANT, BALANCED };

/**
 * @class DeutschJozsa
 * @brief Implements the Deutsch-Jozsa algorithm.
 *
 * Determines with certainty in a single oracle query whether an n-bit boolean
 * function f:{0,1}^n → {0,1} is constant or balanced — exponentially faster
 * than the best classical deterministic algorithm (which needs 2^(n-1)+1 queries).
 *
 * Circuit structure (n+1 qubits, qubit n is the ancilla):
 * @code
 *   |0⟩^⊗n ─── H^⊗n ─── Uf ─── H^⊗n ─── measure q[0..n-1]
 *   |0⟩    ─── X ──── H ─── Uf ───────────────────────────
 * @endcode
 * If all measured qubits are 0 → CONSTANT; any qubit is 1 → BALANCED.
 *
 * @note The oracle @p oracle must implement the phase-kickback form:
 *       |x⟩|−⟩ → (−1)^f(x) |x⟩|−⟩.
 */
class DeutschJozsa {
    private:
        /** @brief Number of input qubits (the ancilla qubit is implicit, at index n). */
        int n;

        /** @brief Internal quantum state of n+1 qubits. */
        QuantumState state;

        /** @brief The quantum oracle implementing f. Called once during run(). */
        std::function<void(QuantumState&)> oracle;

    public:
        /**
         * @brief Constructs a Deutsch-Jozsa instance for an n-bit oracle.
         *
         * @param n      Number of input qubits.
         * @param oracle Quantum oracle implementing f as a phase-kickback unitary.
         */
        DeutschJozsa(int n, std::function<void(QuantumState&)> oracle);

        /**
         * @brief Runs the Deutsch-Jozsa circuit and returns the function classification.
         *
         * @return FunctionType::CONSTANT if all input qubits measure 0,
         *         FunctionType::BALANCED if any input qubit measures 1.
         */
        FunctionType run();
};

#endif /* DEUTSCHJOZSA_HPP_ */
