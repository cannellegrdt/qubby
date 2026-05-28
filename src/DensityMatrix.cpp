/*
 * Project: qubby
 * File name: DensityMatrix.cpp
 * Author: Cannelle Gourdet - lankley
 * File description: Implements DensityMatrix: gate application via left/right multiplication,
 *                   Kraus-operator noise channels (bit-flip, phase-flip, depolarizing),
 *                   and Born-rule measurement on mixed states.
 */

#include "DensityMatrix.hpp"
#include <stdexcept>
#include <cmath>
#include <random>
#include <chrono>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <stdint.h>

/**
 * @details
 * Stores ρ in row-major order: element ρ_{ij} at rho[i * dim + j].
 * The initial pure state |0...0⟩⟨0...0| sets only rho[0] = 1, all others = 0.
 */
void DensityMatrix::initialize(int n) {
    if (n > DM_MAX_QUBITS)
        throw std::runtime_error(
            "Density-matrix noise simulation is limited to " +
            std::to_string(DM_MAX_QUBITS) +
            " qubits (4^n × 16 bytes would exceed 512 MB)");

    num_qubits = n;
    dim = 1 << n;
    rho.assign(static_cast<size_t>(dim) * dim, std::complex<double>(0.0, 0.0));
    rho[0] = 1.0;
}

/**
 * @details
 * The bitwise index trick is the same used in QuantumState::applyGate():
 * for each "half-index" h (ranging over 2^(n-1) values), reconstruct
 * i0 (bit q = 0) and i1 (bit q = 1), then apply the 2×2 gate to
 * every column in one pass.
 */
void DensityMatrix::gateLeft(int q, Matrix U) {
    long long half = dim / 2;
    for (long long h = 0; h < half; h++) {
        uint64_t mask = (1ULL << q) - 1;
        uint64_t i0 = ((h >> q) << (q + 1)) | (h & mask);
        uint64_t i1 = i0 | (1ULL << q);

        for (int j = 0; j < dim; j++) {
            auto a = rho[i0 * dim + j];
            auto b = rho[i1 * dim + j];
            rho[i0 * dim + j] = U.m00 * a + U.m01 * b;
            rho[i1 * dim + j] = U.m10 * a + U.m11 * b;
        }
    }
}

/**
 * @details
 * Right-multiplies by U†: ρ → ρ U_q†.
 * The conjugate-transpose of a 2×2 gate is:
 *   [U†]_{00} = U*_{00},  [U†]_{01} = U*_{10}
 *   [U†]_{10} = U*_{01},  [U†]_{11} = U*_{11}
 * So the update for column pair (j0, j1) is:
 *   ρ'[i,j0] = U*_{00}·ρ[i,j0] + U*_{10}·ρ[i,j1]
 *   ρ'[i,j1] = U*_{01}·ρ[i,j0] + U*_{11}·ρ[i,j1]
 */
void DensityMatrix::gateRight(int q, Matrix U) {
    long long half = dim / 2;
    for (long long h = 0; h < half; h++) {
        uint64_t mask = (1ULL << q) - 1;
        uint64_t j0 = ((h >> q) << (q + 1)) | (h & mask);
        uint64_t j1 = j0 | (1ULL << q);

        for (int i = 0; i < dim; i++) {
            auto a = rho[i * dim + j0];
            auto b = rho[i * dim + j1];
            rho[i * dim + j0] = std::conj(U.m00) * a + std::conj(U.m01) * b;
            rho[i * dim + j1] = std::conj(U.m10) * a + std::conj(U.m11) * b;
        }
    }
}

void DensityMatrix::applyGate(int q, Matrix U) {
    gateLeft(q, U);
    gateRight(q, U);
}

/**
 * @details
 * For each Kraus operator K in ops:
 *   1. Copy ρ into temp.
 *   2. Apply K to temp from the left (same row-pair loop as gateLeft).
 *   3. Accumulate K†·temp into new_rho from the right (same column-pair loop as gateRight).
 *
 * After all operators are processed, set ρ = new_rho.
 * The Kraus completeness condition Σ K†K = I guarantees Tr(ρ') = Tr(ρ) = 1.
 */
void DensityMatrix::kraus(int q, const std::vector<Matrix>& ops) {
    std::vector<std::complex<double>> new_rho(static_cast<size_t>(dim) * dim, 0.0);
    long long half = dim / 2;

    for (const auto& K : ops) {
        auto temp = rho;

        for (long long h = 0; h < half; h++) {
            uint64_t mask = (1ULL << q) - 1;
            uint64_t i0 = ((h >> q) << (q + 1)) | (h & mask);
            uint64_t i1 = i0 | (1ULL << q);

            for (int j = 0; j < dim; j++) {
                auto a = temp[i0 * dim + j];
                auto b = temp[i1 * dim + j];
                temp[i0 * dim + j] = K.m00 * a + K.m01 * b;
                temp[i1 * dim + j] = K.m10 * a + K.m11 * b;
            }
        }

        for (long long h = 0; h < half; h++) {
            uint64_t mask = (1ULL << q) - 1;
            uint64_t j0 = ((h >> q) << (q + 1)) | (h & mask);
            uint64_t j1 = j0 | (1ULL << q);

            for (int i = 0; i < dim; i++) {
                auto a = temp[i * dim + j0];
                auto b = temp[i * dim + j1];
                new_rho[i * dim + j0] += std::conj(K.m00) * a + std::conj(K.m01) * b;
                new_rho[i * dim + j1] += std::conj(K.m10) * a + std::conj(K.m11) * b;
            }
        }
    }

    rho = std::move(new_rho);
}

void DensityMatrix::xGate(int q) {
    Matrix U = {0, 1, 1, 0};
    applyGate(q, U);
}

void DensityMatrix::hGate(int q) {
    double v = 1.0 / std::sqrt(2.0);
    Matrix U = {v, v, v, -v};
    applyGate(q, U);
}

void DensityMatrix::zGate(int q) {
    Matrix U = {1, 0, 0, -1};
    applyGate(q, U);
}

void DensityMatrix::rxGate(double theta, int q) {
    double c = std::cos(theta / 2.0);
    std::complex<double> s = std::complex<double>(0, -1) * std::sin(theta / 2.0);
    Matrix U = {c, s, s, c};
    applyGate(q, U);
}

void DensityMatrix::ryGate(double theta, int q) {
    double c = std::cos(theta / 2.0);
    double s = std::sin(theta / 2.0);
    Matrix U = {c, -s, s, c};
    applyGate(q, U);
}

void DensityMatrix::rzGate(double theta, int q) {
    auto neg = std::exp(std::complex<double>(0, -theta / 2.0));
    auto pos = std::exp(std::complex<double>(0,  theta / 2.0));
    Matrix U = {neg, 0, 0, pos};
    applyGate(q, U);
}

/**
 * @details
 * CNOT left-multiplies row indices: for each pair (i0, i1) differing at bit tgt
 * where the control bit in i0 is 1, swap the entire row pair across all columns.
 * CNOT is its own inverse, so the same operation is used for both left and right.
 */
void DensityMatrix::cnotLeft(int ctrl, int tgt) {
    long long half = dim / 2;
    for (long long h = 0; h < half; h++) {
        uint64_t mask = (1ULL << tgt) - 1;
        uint64_t i0 = ((h >> tgt) << (tgt + 1)) | (h & mask);
        uint64_t i1 = i0 | (1ULL << tgt);

        if ((i0 >> ctrl) & 1) {
            for (int j = 0; j < dim; j++)
                std::swap(rho[i0 * dim + j], rho[i1 * dim + j]);
        }
    }
}

void DensityMatrix::cnotRight(int ctrl, int tgt) {
    long long half = dim / 2;
    for (long long h = 0; h < half; h++) {
        uint64_t mask = (1ULL << tgt) - 1;
        uint64_t j0 = ((h >> tgt) << (tgt + 1)) | (h & mask);
        uint64_t j1 = j0 | (1ULL << tgt);

        if ((j0 >> ctrl) & 1) {
            for (int i = 0; i < dim; i++)
                std::swap(rho[i * dim + j0], rho[i * dim + j1]);
        }
    }
}

void DensityMatrix::cnotGate(int ctrl, int tgt) {
    cnotLeft(ctrl, tgt);
    cnotRight(ctrl, tgt);
}

void DensityMatrix::swapGate(int q0, int q1) {
    cnotGate(q0, q1);
    cnotGate(q1, q0);
    cnotGate(q0, q1);
}

void DensityMatrix::toffoliLeft(int c0, int c1, int tgt) {
    long long half = dim / 2;
    for (long long h = 0; h < half; h++) {
        uint64_t mask = (1ULL << tgt) - 1;
        uint64_t i0 = ((h >> tgt) << (tgt + 1)) | (h & mask);
        uint64_t i1 = i0 | (1ULL << tgt);

        if (((i0 >> c0) & 1) && ((i0 >> c1) & 1)) {
            for (int j = 0; j < dim; j++)
                std::swap(rho[i0 * dim + j], rho[i1 * dim + j]);
        }
    }
}

void DensityMatrix::toffoliRight(int c0, int c1, int tgt) {
    long long half = dim / 2;
    for (long long h = 0; h < half; h++) {
        uint64_t mask = (1ULL << tgt) - 1;
        uint64_t j0 = ((h >> tgt) << (tgt + 1)) | (h & mask);
        uint64_t j1 = j0 | (1ULL << tgt);

        if (((j0 >> c0) & 1) && ((j0 >> c1) & 1)) {
            for (int i = 0; i < dim; i++)
                std::swap(rho[i * dim + j0], rho[i * dim + j1]);
        }
    }
}

void DensityMatrix::toffoliGate(int c0, int c1, int tgt) {
    toffoliLeft(c0, c1, tgt);
    toffoliRight(c0, c1, tgt);
}

/**
 * @details
 * Kraus operators: K0 = √(1-p)·I,  K1 = √p·X.
 * Completeness: K0†K0 + K1†K1 = (1-p)I + pI = I ✓
 */
void DensityMatrix::applyBitFlip(int q, double p) {
    double sq = std::sqrt(1.0 - p);
    double sp = std::sqrt(p);
    std::vector<Matrix> ops = {
        {sq, 0, 0, sq},
        {0, sp, sp, 0}
    };
    kraus(q, ops);
}

/**
 * @details
 * Kraus operators: K0 = √(1-p)·I,  K1 = √p·Z.
 * Completeness: K0†K0 + K1†K1 = (1-p)I + pI = I ✓
 */
void DensityMatrix::applyPhaseFlip(int q, double p) {
    double sq = std::sqrt(1.0 - p);
    double sp = std::sqrt(p);
    std::vector<Matrix> ops = {
        {sq, 0, 0, sq},
        {sp, 0, 0, -sp}
    };
    kraus(q, ops);
}

/**
 * @details
 * Kraus operators: K0 = √(1-p)·I,  K1 = √(p/3)·X,  K2 = √(p/3)·Y,  K3 = √(p/3)·Z.
 * Completeness: K0†K0 + K1†K1 + K2†K2 + K3†K3 = (1-p)I + 3·(p/3)I = I ✓
 *
 * Y = [[0, -i], [i, 0]], so K2 = √(p/3)·Y has entries:
 *   m01 = -i√(p/3),  m10 = i√(p/3).
 */
void DensityMatrix::applyDepolarizing(int q, double p) {
    double k0 = std::sqrt(1.0 - p);
    double k1 = std::sqrt(p / 3.0);
    std::complex<double> iy(0.0, k1);
    std::vector<Matrix> ops = {
        {k0, 0, 0, k0},
        {0, k1, k1, 0},
        {0, -iy, iy, 0},
        {k1, 0, 0, -k1}
    };
    kraus(q, ops);
}

/**
 * @details
 * For a mixed state, the probability of outcome i is the i-th diagonal element
 * of the density matrix: P(i) = ρ[i,i] (always real and non-negative for a
 * valid density matrix).
 *
 * After sampling, the state is projected onto the pure state |result⟩⟨result|:
 * all off-diagonal and other diagonal elements are zeroed.
 */
int DensityMatrix::measure() {
    std::vector<double> probs(dim);
    for (int i = 0; i < dim; i++)
        probs[i] = std::real(rho[static_cast<size_t>(i) * dim + i]);

    static std::mt19937 gen = []() {
        std::random_device rd;
        auto seed = rd() ^ static_cast<uint32_t>(
            std::chrono::high_resolution_clock::now().time_since_epoch().count());
        return std::mt19937(seed);
    }();
    std::discrete_distribution<> dist(probs.begin(), probs.end());
    int result = dist(gen);

    std::fill(rho.begin(), rho.end(), std::complex<double>(0.0, 0.0));
    rho[static_cast<size_t>(result) * dim + result] = 1.0;
    return result;
}

/**
 * @details
 * Iterates over diagonal elements ρ[i,i] and prints each state with probability > 1e-6.
 * Unlike a pure-state printState(), only probability is shown (not the complex amplitude),
 * since off-diagonal coherences represent quantum correlations, not observable probabilities.
 */
double DensityMatrix::getDiagonal(int i) const {
    return std::real(rho[static_cast<size_t>(i) * dim + i]);
}

void DensityMatrix::printState() {
    for (int i = 0; i < dim; i++) {
        double prob = std::real(rho[static_cast<size_t>(i) * dim + i]);
        if (prob > 1e-6) {
            std::string binary;
            for (int b = num_qubits - 1; b >= 0; b--)
                binary += ((i >> b) & 1) ? '1' : '0';
            std::cout << "|" << std::setw(num_qubits) << i << "⟩"
                      << "  " << binary
                      << "  " << std::fixed << std::setprecision(3)
                      << (prob * 100.0) << "%\n";
        }
    }
}
