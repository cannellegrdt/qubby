/*
 * Project: qubby
 * File name: BlochSphere.hpp
 * Author: Cannelle Gourdet - lankley
 * File description: Declares the BlochSphere class for visualising a single-qubit state
 *                   on the Bloch sphere via ASCII terminal art and SVG export.
 */

#ifndef BLOCHSPHERE_HPP_
    #define BLOCHSPHERE_HPP_

    #include "QuantumState.hpp"
    #include <string>

/**
 * @class BlochSphere
 * @brief Renders a single-qubit reduced state as a Bloch sphere diagram.
 *
 * Works on any qubit in an n-qubit state by computing the reduced density matrix
 * for that qubit (tracing out all others). The Bloch vector (x, y, z) satisfies:
 *   - x = Tr(X ρ_q) = 2 Re(ρ_{01})
 *   - y = Tr(Y ρ_q) = −2 Im(ρ_{01})
 *   - z = Tr(Z ρ_q) = P(q=|0⟩) − P(q=|1⟩)
 *
 * For a pure single-qubit state |r| = 1 (point on the sphere surface).
 * Entangled qubits produce mixed reduced states with |r| < 1 (interior of sphere).
 *
 * @note Only supports pure-state (QuantumState) backends — not density matrix.
 */
class BlochSphere {
public:
    /** @brief Bloch vector and its norm. */
    struct BlochVector {
        double x, y, z;
        double norm;  ///< |r| — 1.0 for pure state, <1.0 for mixed
    };

    /**
     * @brief Computes the Bloch vector for qubit @p qubit from @p state.
     *
     * Traces out all qubits except @p qubit to obtain the 2×2 reduced density
     * matrix, then extracts (x, y, z) via Tr(P ρ_q) for each Pauli P.
     *
     * @param state  The n-qubit quantum state.
     * @param qubit  Index of the qubit to project (0-based).
     * @return The Bloch vector (x, y, z) and its norm.
     * @throws std::runtime_error if @p qubit is out of range.
     */
    static BlochVector compute(const QuantumState& state, int qubit);

    /**
     * @brief Saves an SVG file rendering the Bloch sphere with the state vector.
     *
     * Generates a 420×450 SVG with the sphere outline, dashed equatorial ellipse,
     * labelled axes (x, y, z / |0⟩, |1⟩), and an arrow showing the Bloch vector.
     *
     * @param bv       The Bloch vector to render.
     * @param qubit    Qubit index used in the SVG title.
     * @param filename Output SVG file path.
     * @throws std::runtime_error if the file cannot be opened.
     */
    static void saveSVG(const BlochVector& bv, int qubit, const std::string& filename);

    /**
     * @brief Renders an ASCII Bloch sphere diagram to stdout.
     *
     * Projects the sphere onto the xz-plane (front view). The y-component is shown
     * numerically below the diagram. The state vector endpoint is marked with '*' and
     * the line from the centre is drawn with '.'.
     *
     * @param bv    The Bloch vector to render.
     * @param qubit Qubit index shown in the title.
     */
    static void printASCII(const BlochVector& bv, int qubit);
};

#endif /* BLOCHSPHERE_HPP_ */
