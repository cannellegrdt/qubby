/*
 * Project: qubby
 * File name: Shor.hpp
 * Author: Cannelle Gourdet - lankley
 * File description: Declares the Shor class, which implements Shor's quantum factorisation
 *                   algorithm using QFT-based order finding.
 */

#ifndef SHOR_HPP_
    #define SHOR_HPP_
    #include "QuantumState.hpp"

/**
 * @class Shor
 * @brief Implements Shor's quantum integer factorisation algorithm.
 *
 * Factorises a composite integer N in O((log N)³) quantum gate operations,
 * exponentially faster than the best known classical algorithm.
 *
 * Circuit layout (3n qubits, where n = ⌈log₂ N⌉):
 *   - Qubits [0, 2n-1]: period-finding register (input to QFT).
 *   - Qubits [2n, 3n-1]: modular exponentiation register (ancilla).
 *
 * Steps per attempt:
 * 1. H^⊗2n on the first register → uniform superposition over [0, Q).
 * 2. modularOracle: compute |x⟩|0⟩ → |x⟩|a^x mod N⟩.
 * 3. Measure ancilla qubits → collapse first register to a periodic sub-superposition.
 * 4. QFT on first register → peaks at multiples of Q/r.
 * 5. Measure and extract r via continued-fraction expansion.
 * 6. Compute factors: gcd(a^(r/2) ± 1, N).
 *
 * Up to 10 attempts are made with a random base a each time.
 */
class Shor {
    private:
        /** @brief Internal quantum state (3n qubits per attempt). */
        QuantumState state;

    public:
        /**
         * @brief Runs Shor's algorithm and returns a non-trivial factor pair of N.
         *
         * @param N The composite integer to factorise (N > 3, not a prime power).
         * @return `{p, q}` such that p * q == N, or `{-1, -1}` if all attempts failed.
         */
        std::pair<long long, long long> run(int N);

        /**
         * @brief Applies the modular exponentiation oracle |x⟩|y⟩ → |x⟩|y ⊕ aˣ mod N⟩.
         *
         * Iterates over all basis states, computes aˣ mod N classically via fast
         * modular exponentiation, then XORs the result into the ancilla register.
         * Amplitudes below 1e-15 are skipped for efficiency.
         *
         * @param a Base of the modular exponentiation (2 ≤ a < N, gcd(a, N) == 1).
         * @param N Modulus.
         */
        void modularOracle(long long a, long long N);
};

#endif /* SHOR_HPP_ */
