/*
 * Project: qubby
 * File name: test_quantum_state.cpp
 * Author: Cannelle Gourdet - lankley
 * File description: Criterion unit tests for QuantumState and QuantumCircuit.
 */

#include <criterion/criterion.h>
#include "QuantumState.hpp"
#include "QuantumCircuit.hpp"
#include <cmath>
#include <stdexcept>

static constexpr double EPS = 1e-9;

/* ── initialize ────────────────────────────────────────────────────────────── */

Test(initialize, amplitude_zero_is_one) {
    QuantumState s;
    s.initialize(2);

    cr_assert_float_eq(s.getAmplitude(0).real(), 1.0, EPS);
    cr_assert_float_eq(s.getAmplitude(0).imag(), 0.0, EPS);
}

Test(initialize, other_amplitudes_are_zero) {
    QuantumState s;
    s.initialize(2);

    cr_assert_float_eq(s.getAmplitude(1).real(), 0.0, EPS);
    cr_assert_float_eq(s.getAmplitude(2).real(), 0.0, EPS);
    cr_assert_float_eq(s.getAmplitude(3).real(), 0.0, EPS);
}

Test(initialize, size_matches_two_to_the_n) {
    QuantumState s;
    s.initialize(3);

    cr_assert_eq(s.measure(), 0);
}

Test(initialize, throws_when_memory_exceeded) {
    QuantumState s;
    cr_assert_throw(s.initialize(26), std::runtime_error);
}

/* ── xGate ─────────────────────────────────────────────────────────────────── */

Test(x_gate, flips_ground_state_to_excited) {
    QuantumState s;
    s.initialize(1);
    s.xGate(0);

    cr_assert_float_eq(s.getAmplitude(0).real(), 0.0, EPS);
    cr_assert_float_eq(s.getAmplitude(1).real(), 1.0, EPS);
}

Test(x_gate, applied_twice_is_identity) {
    QuantumState s;
    s.initialize(1);
    s.xGate(0);
    s.xGate(0);

    cr_assert_float_eq(s.getAmplitude(0).real(), 1.0, EPS);
    cr_assert_float_eq(s.getAmplitude(1).real(), 0.0, EPS);
}

Test(x_gate, targets_correct_qubit_in_two_qubit_system) {
    QuantumState s;
    s.initialize(2);
    s.xGate(1);

    cr_assert_float_eq(s.getAmplitude(0).real(), 0.0, EPS);
    cr_assert_float_eq(s.getAmplitude(2).real(), 1.0, EPS);
}

Test(x_gate, deterministic_measure_returns_one) {
    QuantumState s;
    s.initialize(1);
    s.xGate(0);

    cr_assert_eq(s.measure(), 1);
}

/* ── hGate ─────────────────────────────────────────────────────────────────── */

Test(h_gate, creates_equal_superposition) {
    QuantumState s;
    s.initialize(1);
    s.hGate(0);

    double expected = 1.0 / std::sqrt(2.0);
    cr_assert_float_eq(s.getAmplitude(0).real(), expected, EPS);
    cr_assert_float_eq(s.getAmplitude(1).real(), expected, EPS);
}

Test(h_gate, applied_twice_is_identity) {
    QuantumState s;
    s.initialize(1);
    s.hGate(0);
    s.hGate(0);

    cr_assert_float_eq(s.getAmplitude(0).real(), 1.0, EPS);
    cr_assert_float_eq(s.getAmplitude(1).real(), 0.0, EPS);
}

Test(h_gate, h_squared_measure_is_ground_state) {
    QuantumState s;
    s.initialize(1);
    s.hGate(0);
    s.hGate(0);

    cr_assert_eq(s.measure(), 0);
}

/* ── zGate ─────────────────────────────────────────────────────────────────── */

Test(z_gate, leaves_ground_state_unchanged) {
    QuantumState s;
    s.initialize(1);
    s.zGate(0);

    cr_assert_float_eq(s.getAmplitude(0).real(),  1.0, EPS);
    cr_assert_float_eq(s.getAmplitude(1).real(),  0.0, EPS);
}

Test(z_gate, negates_excited_state_amplitude) {
    QuantumState s;
    s.initialize(1);
    s.xGate(0);
    s.zGate(0);

    cr_assert_float_eq(s.getAmplitude(0).real(),  0.0, EPS);
    cr_assert_float_eq(s.getAmplitude(1).real(), -1.0, EPS);
}

Test(z_gate, hzh_equals_x) {
    QuantumState s;
    s.initialize(1);
    s.hGate(0);
    s.zGate(0);
    s.hGate(0);

    cr_assert_float_eq(s.getAmplitude(0).real(), 0.0, EPS);
    cr_assert_float_eq(s.getAmplitude(1).real(), 1.0, EPS);
}

/* ── rxGate ─────────────────────────────────────────────────────────────────── */

Test(rx_gate, zero_angle_is_identity) {
    QuantumState s;
    s.initialize(1);
    s.rxGate(0.0, 0);

    cr_assert_float_eq(s.getAmplitude(0).real(), 1.0, EPS);
    cr_assert_float_eq(s.getAmplitude(0).imag(), 0.0, EPS);
    cr_assert_float_eq(s.getAmplitude(1).real(), 0.0, EPS);
    cr_assert_float_eq(s.getAmplitude(1).imag(), 0.0, EPS);
}

Test(rx_gate, pi_rotation_flips_ground_state) {
    QuantumState s;
    s.initialize(1);
    s.rxGate(M_PI, 0);

    cr_assert_float_eq(s.getAmplitude(0).real(),  0.0, EPS);
    cr_assert_float_eq(s.getAmplitude(0).imag(),  0.0, EPS);
    cr_assert_float_eq(s.getAmplitude(1).real(),  0.0, EPS);
    cr_assert_float_eq(s.getAmplitude(1).imag(), -1.0, EPS);
}

Test(rx_gate, pi_rotation_measures_one) {
    QuantumState s;
    s.initialize(1);
    s.rxGate(M_PI, 0);

    cr_assert_eq(s.measure(), 1);
}

Test(rx_gate, half_pi_creates_superposition) {
    QuantumState s;
    s.initialize(1);
    s.rxGate(M_PI / 2.0, 0);

    double expected = 1.0 / std::sqrt(2.0);
    cr_assert_float_eq(s.getAmplitude(0).real(),  expected, EPS);
    cr_assert_float_eq(s.getAmplitude(0).imag(),  0.0,      EPS);
    cr_assert_float_eq(s.getAmplitude(1).real(),  0.0,      EPS);
    cr_assert_float_eq(s.getAmplitude(1).imag(), -expected, EPS);
}

Test(rx_gate, two_pi_rotations_return_minus_identity) {
    QuantumState s;
    s.initialize(1);
    s.rxGate(M_PI, 0);
    s.rxGate(M_PI, 0);

    cr_assert_float_eq(s.getAmplitude(0).real(), -1.0, EPS);
    cr_assert_float_eq(s.getAmplitude(0).imag(),  0.0, EPS);
    cr_assert_float_eq(s.getAmplitude(1).real(),  0.0, EPS);
    cr_assert_float_eq(s.getAmplitude(1).imag(),  0.0, EPS);
}

/* ── ryGate ─────────────────────────────────────────────────────────────────── */

Test(ry_gate, zero_angle_is_identity) {
    QuantumState s;
    s.initialize(1);
    s.ryGate(0.0, 0);

    cr_assert_float_eq(s.getAmplitude(0).real(), 1.0, EPS);
    cr_assert_float_eq(s.getAmplitude(0).imag(), 0.0, EPS);
    cr_assert_float_eq(s.getAmplitude(1).real(), 0.0, EPS);
    cr_assert_float_eq(s.getAmplitude(1).imag(), 0.0, EPS);
}

Test(ry_gate, pi_rotation_flips_ground_state) {
    QuantumState s;
    s.initialize(1);
    s.ryGate(M_PI, 0);

    cr_assert_float_eq(s.getAmplitude(0).real(), 0.0, EPS);
    cr_assert_float_eq(s.getAmplitude(0).imag(), 0.0, EPS);
    cr_assert_float_eq(s.getAmplitude(1).real(), 1.0, EPS);
    cr_assert_float_eq(s.getAmplitude(1).imag(), 0.0, EPS);
}

Test(ry_gate, pi_rotation_measures_one) {
    QuantumState s;
    s.initialize(1);
    s.ryGate(M_PI, 0);

    cr_assert_eq(s.measure(), 1);
}

Test(ry_gate, half_pi_creates_real_superposition) {
    QuantumState s;
    s.initialize(1);
    s.ryGate(M_PI / 2.0, 0);

    double expected = 1.0 / std::sqrt(2.0);
    cr_assert_float_eq(s.getAmplitude(0).real(),  expected, EPS);
    cr_assert_float_eq(s.getAmplitude(0).imag(),  0.0,      EPS);
    cr_assert_float_eq(s.getAmplitude(1).real(),  expected, EPS);
    cr_assert_float_eq(s.getAmplitude(1).imag(),  0.0,      EPS);
}

Test(ry_gate, two_pi_rotations_return_minus_identity) {
    QuantumState s;
    s.initialize(1);
    s.ryGate(M_PI, 0);
    s.ryGate(M_PI, 0);

    cr_assert_float_eq(s.getAmplitude(0).real(), -1.0, EPS);
    cr_assert_float_eq(s.getAmplitude(0).imag(),  0.0, EPS);
    cr_assert_float_eq(s.getAmplitude(1).real(),  0.0, EPS);
    cr_assert_float_eq(s.getAmplitude(1).imag(),  0.0, EPS);
}

/* ── rzGate ─────────────────────────────────────────────────────────────────── */

Test(rz_gate, zero_angle_is_identity) {
    QuantumState s;
    s.initialize(1);
    s.rzGate(0.0, 0);

    cr_assert_float_eq(s.getAmplitude(0).real(), 1.0, EPS);
    cr_assert_float_eq(s.getAmplitude(0).imag(), 0.0, EPS);
    cr_assert_float_eq(s.getAmplitude(1).real(), 0.0, EPS);
    cr_assert_float_eq(s.getAmplitude(1).imag(), 0.0, EPS);
}

Test(rz_gate, leaves_ground_state_probability_unchanged) {
    QuantumState s;
    s.initialize(1);
    s.rzGate(M_PI / 2.0, 0);

    cr_assert_eq(s.measure(), 0);
}

Test(rz_gate, pi_rotation_adds_phase_to_ground_state) {
    QuantumState s;
    s.initialize(1);
    s.rzGate(M_PI, 0);

    cr_assert_float_eq(s.getAmplitude(0).real(),  0.0, EPS);
    cr_assert_float_eq(s.getAmplitude(0).imag(), -1.0, EPS);
    cr_assert_float_eq(s.getAmplitude(1).real(),  0.0, EPS);
    cr_assert_float_eq(s.getAmplitude(1).imag(),  0.0, EPS);
}

Test(rz_gate, pi_rotation_adds_phase_to_excited_state) {
    QuantumState s;
    s.initialize(1);
    s.xGate(0);
    s.rzGate(M_PI, 0);

    cr_assert_float_eq(s.getAmplitude(0).real(), 0.0, EPS);
    cr_assert_float_eq(s.getAmplitude(0).imag(), 0.0, EPS);
    cr_assert_float_eq(s.getAmplitude(1).real(), 0.0, EPS);
    cr_assert_float_eq(s.getAmplitude(1).imag(), 1.0, EPS);
}

Test(rz_gate, hrzh_equals_x_up_to_phase) {
    QuantumState s;
    s.initialize(1);
    s.hGate(0);
    s.rzGate(M_PI, 0);
    s.hGate(0);

    cr_assert_eq(s.measure(), 1);
}

/* ── cnotGate ──────────────────────────────────────────────────────────────── */

Test(cnot_gate, control_zero_leaves_target_unchanged) {
    QuantumState s;
    s.initialize(2);
    s.cnotGate(0, 1);

    cr_assert_float_eq(s.getAmplitude(0).real(), 1.0, EPS);
    cr_assert_float_eq(s.getAmplitude(1).real(), 0.0, EPS);
    cr_assert_float_eq(s.getAmplitude(2).real(), 0.0, EPS);
    cr_assert_float_eq(s.getAmplitude(3).real(), 0.0, EPS);
}

Test(cnot_gate, control_one_flips_target) {
    QuantumState s;
    s.initialize(2);
    s.xGate(0);
    s.cnotGate(0, 1);

    cr_assert_float_eq(s.getAmplitude(1).real(), 0.0, EPS);
    cr_assert_float_eq(s.getAmplitude(3).real(), 1.0, EPS);
}

Test(cnot_gate, applied_twice_is_identity) {
    QuantumState s;
    s.initialize(2);
    s.xGate(0);
    s.cnotGate(0, 1);
    s.cnotGate(0, 1);

    cr_assert_float_eq(s.getAmplitude(1).real(), 1.0, EPS);
    cr_assert_float_eq(s.getAmplitude(3).real(), 0.0, EPS);
}

Test(cnot_gate, creates_bell_state) {
    QuantumState s;
    s.initialize(2);
    s.hGate(0);
    s.cnotGate(0, 1);

    double expected = 1.0 / std::sqrt(2.0);
    cr_assert_float_eq(s.getAmplitude(0).real(), expected, EPS);
    cr_assert_float_eq(s.getAmplitude(3).real(), expected, EPS);
    cr_assert_float_eq(s.getAmplitude(1).real(), 0.0, EPS);
    cr_assert_float_eq(s.getAmplitude(2).real(), 0.0, EPS);
}

/* ── swapGate ──────────────────────────────────────────────────────────────── */

Test(swap_gate, leaves_ground_state_unchanged) {
    QuantumState s;
    s.initialize(2);
    s.swapGate(0, 1);

    cr_assert_eq(s.measure(), 0);
}

Test(swap_gate, swaps_qubit_zero_to_qubit_one) {
    QuantumState s;
    s.initialize(2);
    s.xGate(0);
    s.swapGate(0, 1);

    cr_assert_float_eq(s.getAmplitude(2).real(), 1.0, EPS);
    cr_assert_eq(s.measure(), 2);
}

Test(swap_gate, swaps_qubit_one_to_qubit_zero) {
    QuantumState s;
    s.initialize(2);
    s.xGate(1);
    s.swapGate(0, 1);

    cr_assert_float_eq(s.getAmplitude(1).real(), 1.0, EPS);
    cr_assert_eq(s.measure(), 1);
}

Test(swap_gate, applied_twice_is_identity) {
    QuantumState s;
    s.initialize(2);
    s.xGate(0);
    s.swapGate(0, 1);
    s.swapGate(0, 1);

    cr_assert_float_eq(s.getAmplitude(1).real(), 1.0, EPS);
    cr_assert_eq(s.measure(), 1);
}

/* ── toffoliGate ───────────────────────────────────────────────────────────── */

Test(toffoli_gate, both_controls_zero_leaves_target_unchanged) {
    QuantumState s;
    s.initialize(3);
    s.toffoliGate(0, 1, 2);

    cr_assert_eq(s.measure(), 0);
}

Test(toffoli_gate, single_control_set_leaves_target_unchanged) {
    QuantumState s;
    s.initialize(3);
    s.xGate(0);
    s.toffoliGate(0, 1, 2);

    cr_assert_eq(s.measure(), 1);
}

Test(toffoli_gate, both_controls_set_flips_target) {
    QuantumState s;
    s.initialize(3);
    s.xGate(0);
    s.xGate(1);
    s.toffoliGate(0, 1, 2);

    cr_assert_float_eq(s.getAmplitude(7).real(), 1.0, EPS);
    cr_assert_eq(s.measure(), 7);
}

Test(toffoli_gate, applied_twice_is_identity) {
    QuantumState s;
    s.initialize(3);
    s.xGate(0);
    s.xGate(1);
    s.toffoliGate(0, 1, 2);
    s.toffoliGate(0, 1, 2);

    cr_assert_eq(s.measure(), 3);
}

/* ── measure ───────────────────────────────────────────────────────────────── */

Test(measure, ground_state_always_returns_zero) {
    QuantumState s;
    s.initialize(3);

    cr_assert_eq(s.measure(), 0);
}

Test(measure, excited_state_always_returns_one) {
    QuantumState s;
    s.initialize(1);
    s.xGate(0);

    cr_assert_eq(s.measure(), 1);
}

Test(measure, collapses_superposition_to_same_result) {
    QuantumState s;
    s.initialize(1);
    s.hGate(0);

    int first  = s.measure();
    int second = s.measure();
    cr_assert_eq(first, second);
}

Test(measure, bell_state_result_is_zero_or_three) {
    QuantumState s;
    s.initialize(2);
    s.hGate(0);
    s.cnotGate(0, 1);

    int result = s.measure();
    cr_assert(result == 0 || result == 3,
              "Bell state measurement must be 0 or 3, got %d", result);
}

/* ── QuantumCircuit ────────────────────────────────────────────────────────── */

Test(quantum_circuit, throws_on_missing_file) {
    cr_assert_throw(
        QuantumCircuit(1).load("nonexistent_file.qasm"),
        std::runtime_error
    );
}

Test(quantum_circuit, throws_on_unknown_gate) {
    cr_assert_throw(
        QuantumCircuit(1).load("tests/functional/circuits/invalid_gate.qasm"),
        std::runtime_error
    );
}

Test(quantum_circuit, x_gate_circuit_measures_one) {
    QuantumCircuit c(1);
    c.load("tests/functional/circuits/x_gate.qasm");

    cr_assert_eq(c.measure(), 1);
}

Test(quantum_circuit, hh_identity_circuit_measures_zero) {
    QuantumCircuit c(1);
    c.load("tests/functional/circuits/hh_identity.qasm");

    cr_assert_eq(c.measure(), 0);
}

Test(quantum_circuit, comments_are_ignored) {
    QuantumCircuit c(1);
    c.load("tests/functional/circuits/comments.qasm");

    cr_assert_eq(c.measure(), 1);
}

Test(quantum_circuit, z_gate_circuit_on_ground_state_measures_zero) {
    QuantumCircuit c(1);
    c.load("tests/functional/circuits/z_gate.qasm");

    cr_assert_eq(c.measure(), 0);
}

Test(quantum_circuit, throws_on_malformed_qubit_token) {
    cr_assert_throw(
        QuantumCircuit(1).load("tests/functional/circuits/malformed_qubit.qasm"),
        std::runtime_error
    );
}

Test(quantum_circuit, bell_state_circuit_result_is_zero_or_three) {
    QuantumCircuit c(2);
    c.load("tests/functional/circuits/bell_state.qasm");

    int result = c.measure();
    cr_assert(result == 0 || result == 3,
              "Bell state measurement must be 0 or 3, got %d", result);
}

Test(quantum_circuit, rx_pi_circuit_measures_one) {
    QuantumCircuit c(1);
    c.load("tests/functional/circuits/rx_pi.qasm");

    cr_assert_eq(c.measure(), 1);
}

Test(quantum_circuit, ry_pi_circuit_measures_one) {
    QuantumCircuit c(1);
    c.load("tests/functional/circuits/ry_pi.qasm");

    cr_assert_eq(c.measure(), 1);
}

Test(quantum_circuit, rz_circuit_on_ground_state_measures_zero) {
    QuantumCircuit c(1);
    c.load("tests/functional/circuits/rz_ground.qasm");

    cr_assert_eq(c.measure(), 0);
}

Test(quantum_circuit, block_comment_is_ignored) {
    QuantumCircuit c(1);
    c.load("tests/functional/circuits/block_comment.qasm");

    cr_assert_eq(c.measure(), 1);
}

Test(quantum_circuit, qreg_reinitializes_state) {
    QuantumCircuit c(1);
    c.load("tests/functional/circuits/qreg.qasm");

    cr_assert_eq(c.measure(), 1);
}

Test(quantum_circuit, symbolic_pi_angle_measures_one) {
    QuantumCircuit c(1);
    c.load("tests/functional/circuits/rx_symbolic_pi.qasm");

    cr_assert_eq(c.measure(), 1);
}

Test(quantum_circuit, symbolic_pi_div_2_creates_superposition) {
    QuantumCircuit c(1);
    c.load("tests/functional/circuits/ry_symbolic_pi_div2.qasm");

    int result = c.measure();
    cr_assert(result == 0 || result == 1,
              "Ry(pi/2) measure must be 0 or 1, got %d", result);
}

Test(quantum_circuit, multiline_block_comment_is_ignored) {
    QuantumCircuit c(1);
    c.load("tests/functional/circuits/multiline_block_comment.qasm");

    cr_assert_eq(c.measure(), 1);
}

Test(quantum_circuit, symbolic_pi_mult_angle_measures_zero) {
    QuantumCircuit c(1);
    c.load("tests/functional/circuits/rx_symbolic_pi_mult.qasm");

    cr_assert_eq(c.measure(), 0);
}

Test(quantum_circuit, throws_on_malformed_parametric_angle) {
    cr_assert_throw(
        QuantumCircuit(1).load("tests/functional/circuits/malformed_rx_angle.qasm"),
        std::runtime_error
    );
}
