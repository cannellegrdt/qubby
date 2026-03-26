/*
 * Project: qubby
 * File name: VariationalQuantumEigensolver.cpp
 * Author: Cannelle Gourdet - lankley
 * File description: Implements the VQE optimisation loop: ansatz preparation, Pauli expectation
 *                   estimation, and parameter-shift gradient descent.
 */

#include "VariationalQuantumEigensolver.hpp"

VQE::VQE(int n_qubits, std::vector<PauliTerm> hamiltonian)
    : n_qubits(n_qubits), hamiltonian(std::move(hamiltonian)) {}

/**
 * @details
 * Hardware-efficient ansatz: one layer of Ry rotations followed by a linear
 * CNOT entangling chain. This produces a 1-design that can represent a broad
 * family of ground states while keeping the circuit depth linear in n.
 */
void VQE::prepareAnsatz(QuantumState &state, const std::vector<double> &theta) const {
    state.initialize(n_qubits);
    for (int i=0; i<n_qubits; i++) state.ryGate(theta[i], i);
    for (int i=0; i<n_qubits-1; i++) state.cnotGate(i, i+1);
}

/**
 * @details
 * Measurement basis rotation before sampling:
 *   - PauliOp::X: apply H to rotate X eigenbasis → Z eigenbasis.
 *   - PauliOp::Y: apply Rx(π/2) to rotate Y eigenbasis → Z eigenbasis.
 *   - PauliOp::Z / I: no rotation needed.
 *
 * The expectation is then computed analytically from amplitudes:
 *   ⟨P⟩ = Σₓ |ψ(x)|² · (−1)^(number of active non-I bits set in x)
 * where "active" means the qubit's operator is not I.
 * The result is scaled by term.coeff before returning.
 */
double VQE::measurePauliTerm(const PauliTerm &term, const std::vector<double> &theta) const {
    QuantumState state;
    prepareAnsatz(state, theta);

    for (int i=0; i<n_qubits; i++) {
        if (term.ops[i] == PauliOp::X) state.hGate(i);
        else if (term.ops[i] == PauliOp::Y) state.rxGate(M_PI/2.0, i);
    }

    int dim = 1 << n_qubits;
    double expectation = 0.0;
    for (int x=0; x<dim; x++) {
        double prob = std::norm(state.getAmplitude(x));
        int sign = 1;
        for (int i=0; i<n_qubits; i++) {
            if (term.ops[i] != PauliOp::I && ((x >> i) & 1)) sign *= -1;
        }
        expectation += prob * sign;
    }
    return term.coeff * expectation;
}

double VQE::expectationValue(const std::vector<double> &theta) const {
    double energy = 0.0;
    for (const auto &term : hamiltonian)
        energy += measurePauliTerm(term, theta);
    return energy;
}

std::vector<double> VQE::getOptimalParams() const {
    return params;
}

/**
 * @details
 * Optimisation loop (up to 500 iterations, Adam optimiser):
 * 1. For each parameter θᵢ, evaluate E(θᵢ + π/2) and E(θᵢ - π/2) to get
 *    ∂E/∂θᵢ via the parameter-shift rule (exact gradient, not finite difference).
 * 2. Compute gradient L2 norm; break early if < 1e-6.
 * 3. Adam update with β₁=0.9, β₂=0.999, ε=1e-8, lr=0.01:
 *      mᵢ ← β₁·mᵢ + (1−β₁)·gᵢ
 *      vᵢ ← β₂·vᵢ + (1−β₂)·gᵢ²
 *      θᵢ ← θᵢ − lr · m̂ᵢ / (√v̂ᵢ + ε)   (bias-corrected moments)
 *
 * The parameter-shift rule is exact for gates of the form e^(-iθP/2) (e.g. Ry),
 * so the gradient is exact, not an approximation.
 */
double VQE::run() {
    int n_params = n_qubits;
    params.assign(n_params, 0.1);

    const double lr    = 0.01;
    const double beta1 = 0.9;
    const double beta2 = 0.999;
    const double eps   = 1e-8;
    const double grad_tol = 1e-6;

    std::vector<double> m(n_params, 0.0);
    std::vector<double> v(n_params, 0.0);

    for (int iter = 1; iter <= 500; iter++) {
        std::vector<double> grad(n_params);
        for (int i = 0; i < n_params; i++) {
            auto p_plus = params, p_minus = params;
            p_plus[i]  += M_PI / 2.0;
            p_minus[i] -= M_PI / 2.0;
            grad[i] = (expectationValue(p_plus) - expectationValue(p_minus)) / 2.0;
        }

        double grad_norm = 0.0;
        for (int i = 0; i < n_params; i++)
            grad_norm += grad[i] * grad[i];
        if (std::sqrt(grad_norm) < grad_tol)
            break;

        for (int i = 0; i < n_params; i++) {
            m[i] = beta1 * m[i] + (1.0 - beta1) * grad[i];
            v[i] = beta2 * v[i] + (1.0 - beta2) * grad[i] * grad[i];

            double m_hat = m[i] / (1.0 - std::pow(beta1, iter));
            double v_hat = v[i] / (1.0 - std::pow(beta2, iter));

            params[i] -= lr * m_hat / (std::sqrt(v_hat) + eps);
        }
    }
    return expectationValue(params);
}
