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
