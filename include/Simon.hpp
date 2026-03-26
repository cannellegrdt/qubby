/*
 * Project: qubby
 * File name: Simon.hpp
 * Author: Cannelle Gourdet - lankley
 * File description: Declares the Simon class, which implements Simon's algorithm to find
 *                   the hidden period of a 2-to-1 boolean function in a single oracle pass.
 */

#ifndef SIMON_HPP_
    #define SIMON_HPP_
    #include "QuantumState.hpp"
    #include <functional>

/**
 * @class Simon
 * @brief Implements Simon's algorithm to find the hidden period of a 2-to-1 function.
 *
 * Given f:{0,1}^n → {0,1}^n satisfying f(x) = f(x ⊕ s) for some secret string s,
 * Simon's algorithm finds s using O(n) oracle calls — exponentially faster than
 * any classical algorithm (which needs Ω(2^(n/2)) calls).
 *
 * Each iteration produces a uniformly random y ∈ {0,1}^n satisfying y·s = 0 (mod 2).
 * After n-1 linearly independent equations are collected, s is recovered by GF(2) Gaussian
 * elimination.
 *
 * Circuit per iteration (2n qubits, upper n = input, lower n = output):
 * @code
 *   |0⟩^⊗n ─── H^⊗n ─── Uf ─── measure lower ─── H^⊗n ─── measure upper
 *   |0⟩^⊗n ──────────── Uf ─── measure lower
 * @endcode
 */
class Simon {
    private:
        /** @brief Internal quantum state (2n qubits per iteration). */
        QuantumState state;

        /**
         * @brief Recovers s from a system of GF(2) linear equations y·s = 0.
         *
         * Performs Gaussian elimination over GF(2) on the collected measurement
         * vectors. Free variables are set to 1 to produce a non-trivial solution.
         *
         * @param equations List of measurement bit-vectors (each as an int bitmask).
         * @param n         Number of bits (dimension of the GF(2) space).
         * @return The secret string s as an integer bitmask.
         */
        int solveGF2(std::vector<int> &equations, int n);

    public:
        /**
         * @brief Runs Simon's algorithm and returns the hidden period s.
         *
         * Collects n-1 linearly independent measurement outcomes via repeated
         * quantum circuit executions, then calls solveGF2() to extract s.
         *
         * @param n Number of input bits.
         * @param f The 2-to-1 oracle function f:{0,1}^n → {0,1}^n.
         * @return The hidden period s as an integer bitmask.
         */
        int run(int n, std::function<int(int)> f);

        /**
         * @brief Applies the Simon oracle |x⟩|y⟩ → |x⟩|y ⊕ f(x)⟩.
         *
         * Iterates over all 2^(2n) basis states, extracts the upper n bits as x,
         * computes f(x), and XORs it into the lower n bits. Amplitudes below 1e-12
         * are skipped.
         *
         * @param f Oracle function mapping x → f(x) (integer bitmasks).
         * @param n Number of input/output bits.
         */
        void applyOracle(std::function<int(int)> f, int n);
};

#endif /* SIMON_HPP_ */
