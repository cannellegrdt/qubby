/*
 * Project: qubby
 * File name: Matrix.hpp
 * Author: Cannelle Gourdet - lankley
 * File description: ...
 */

#ifndef MATRIX_HPP_
    #define MATRIX_HPP_

    #include <complex>

/**
 * @struct Matrix
 * @brief Represents a 2×2 complex matrix used as a single-qubit gate.
 *
 * A unitary 2×2 matrix acting on the |0⟩/|1⟩ subspace of one qubit:
 * @code
 *   | m00  m01 |
 *   | m10  m11 |
 * @endcode
 */
struct Matrix {
    std::complex<double> m00, m01;
    std::complex<double> m10, m11;
};

#endif /* MATRIX_HPP_ */
