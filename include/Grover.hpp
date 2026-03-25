/*
 * Project: qubby
 * File name: Grover.hpp
 * Author: Cannelle Gourdet - lankley
 * File description: Declares the Grover class, which implements Grover's amplitude-amplification
 *                   algorithm to search an unstructured database in O(√N) oracle queries.
 */

#ifndef GROVER_HPP_
    #define GROVER_HPP_
    #include "QuantumState.hpp"

/**
 * @class Grover
 * @brief Implements Grover's search algorithm.
 *
 * Searches for a marked basis state @p target in a 2^n-element unstructured search space
 * using O(√(2^n)) iterations — quadratically faster than any classical algorithm.
 *
 * Each iteration applies two operators:
 * 1. **Oracle** (phaseFlip): marks the target state with a -1 phase.
 * 2. **Diffusion** (inversion about the mean): amplifies the target's amplitude.
 *
 * After k ≈ (π/4)√(2^n) iterations, a measurement returns @p target with high probability.
 */
class Grover {
    private:
        /** @brief Number of qubits (search space has 2^n_qubits elements). */
        int n_qubits;

        /** @brief Index of the basis state to search for. */
        int target;

        /** @brief Internal quantum state used throughout the algorithm. */
        QuantumState state;

        /**
         * @brief Marks the target state by flipping its amplitude sign (phase kickback).
         *
         * Calls QuantumState::phaseFlip(target): ψ[target] *= -1.
         */
        void applyOracle();

        /**
         * @brief Applies the Grover diffusion operator (inversion about the mean).
         *
         * Sequence: H^⊗n → phaseFlipAllExceptZero → H^⊗n.
         * Geometrically reflects all amplitudes about their mean, amplifying the
         * marked state's amplitude after each oracle call.
         */
        void applyDiffusion();

    public:
        /**
         * @brief Constructs a Grover search instance.
         *
         * @param n_qubits Number of qubits (search space size = 2^n_qubits).
         * @param target   Basis-state index of the element to search for.
         */
        Grover(int n_qubits, int target);

        /**
         * @brief Runs the Grover algorithm and returns the measured result.
         *
         * Initialises the state to |00...0⟩, applies H^⊗n, then iterates
         * oracle + diffusion for k ≈ (π/4)√(2^n) steps, and finally measures.
         *
         * @return The measured basis-state index — @p target with high probability.
         */
        int run();
};

#endif /* GROVER_HPP_ */
