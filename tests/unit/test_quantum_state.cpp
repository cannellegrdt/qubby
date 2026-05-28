/*
 * Project: qubby
 * File name: test_quantum_state.cpp
 * Author: Cannelle Gourdet - lankley
 * File description: Criterion unit tests for QuantumState and QuantumCircuit.
 */

#include <criterion/criterion.h>
#include <criterion/redirect.h>
#include "QuantumState.hpp"
#include "QuantumCircuit.hpp"
#include "DensityMatrix.hpp"
#include "QuantumFourierTransform.hpp"
#include "Grover.hpp"
#include "DeutschJozsa.hpp"
#include "Shor.hpp"
#include "Simon.hpp"
#include "VariationalQuantumEigensolver.hpp"
#include <cmath>
#include <stdexcept>
#include <numeric>

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

/* ── measureQubit ──────────────────────────────────────────────────────────── */

Test(measure_qubit, ground_state_returns_zero) {
    QuantumState s;
    s.initialize(1);

    cr_assert_eq(s.measureQubit(0), 0);
}

Test(measure_qubit, excited_state_returns_one) {
    QuantumState s;
    s.initialize(1);
    s.xGate(0);

    cr_assert_eq(s.measureQubit(0), 1);
}

Test(measure_qubit, collapses_the_qubit) {
    QuantumState s;
    s.initialize(1);
    s.hGate(0);

    int first  = s.measureQubit(0);
    int second = s.measureQubit(0);
    cr_assert_eq(first, second);
}

Test(measure_qubit, only_target_qubit_is_affected) {
    QuantumState s;
    s.initialize(2);
    s.xGate(1);

    cr_assert_eq(s.measureQubit(0), 0);
    cr_assert_eq(s.measureQubit(1), 1);
}

/* ── phaseFlip ─────────────────────────────────────────────────────────────── */

Test(phase_flip, negates_amplitude_at_index) {
    QuantumState s;
    s.initialize(1);
    s.phaseFlip(0);

    cr_assert_float_eq(s.getAmplitude(0).real(), -1.0, EPS);
    cr_assert_float_eq(s.getAmplitude(0).imag(),  0.0, EPS);
}

Test(phase_flip, applied_twice_is_identity) {
    QuantumState s;
    s.initialize(1);
    s.phaseFlip(0);
    s.phaseFlip(0);

    cr_assert_float_eq(s.getAmplitude(0).real(), 1.0, EPS);
}

Test(phase_flip, does_not_change_measurement_probability) {
    QuantumState s;
    s.initialize(1);
    s.phaseFlip(0);

    cr_assert_eq(s.measure(), 0);
}

/* ── phaseFlipAllExceptZero ────────────────────────────────────────────────── */

Test(phase_flip_all_except_zero, ground_state_is_preserved) {
    QuantumState s;
    s.initialize(1);
    s.phaseFlipAllExceptZero();

    cr_assert_float_eq(s.getAmplitude(0).real(), 1.0, EPS);
    cr_assert_float_eq(s.getAmplitude(1).real(), 0.0, EPS);
}

Test(phase_flip_all_except_zero, negates_non_zero_index_amplitudes) {
    QuantumState s;
    s.initialize(1);
    s.hGate(0);
    s.phaseFlipAllExceptZero();

    double expected = 1.0 / std::sqrt(2.0);
    cr_assert_float_eq(s.getAmplitude(0).real(),  expected, EPS);
    cr_assert_float_eq(s.getAmplitude(1).real(), -expected, EPS);
}

/* ── applyControlledRz ─────────────────────────────────────────────────────── */

Test(controlled_rz, no_effect_when_control_is_zero) {
    QuantumState s;
    s.initialize(2);
    s.applyControlledRz(0, 1, M_PI);

    cr_assert_float_eq(s.getAmplitude(0).real(), 1.0, EPS);
    cr_assert_float_eq(s.getAmplitude(0).imag(), 0.0, EPS);
}

Test(controlled_rz, applies_phase_when_control_is_one) {
    QuantumState s;
    s.initialize(2);
    s.xGate(0);
    s.xGate(1);
    s.applyControlledRz(0, 1, M_PI);

    cr_assert_float_eq(s.getAmplitude(3).real(), -1.0, EPS);
    cr_assert_float_eq(s.getAmplitude(3).imag(),  0.0, EPS);
}

Test(controlled_rz, applies_quarter_turn_phase) {
    QuantumState s;
    s.initialize(2);
    s.xGate(0);
    s.xGate(1);
    s.applyControlledRz(0, 1, M_PI / 2.0);

    cr_assert_float_eq(s.getAmplitude(3).real(), 0.0, EPS);
    cr_assert_float_eq(s.getAmplitude(3).imag(), 1.0, EPS);
}

/* ── QFT ───────────────────────────────────────────────────────────────────── */

Test(qft, one_qubit_equals_hadamard_on_ground_state) {
    QFT qft(1);
    qft.run();

    QuantumState& state = qft.getState();
    double expected = 1.0 / std::sqrt(2.0);
    cr_assert_float_eq(state.getAmplitude(0).real(), expected, EPS);
    cr_assert_float_eq(state.getAmplitude(1).real(), expected, EPS);
}

Test(qft, two_qubits_produces_uniform_superposition_on_ground_state) {
    QFT qft(2);
    qft.run();

    QuantumState& state = qft.getState();
    double expected = 0.5;
    for (int i = 0; i < 4; i++)
        cr_assert_float_eq(std::abs(state.getAmplitude(i)), expected, EPS);
}

/* ── Grover ────────────────────────────────────────────────────────────────── */

Test(grover, returns_valid_index_in_range) {
    int result = Grover(2, 1).run();
    cr_assert(result >= 0 && result < 4,
              "Grover(2,1) must return an index in [0,4), got %d", result);
}

Test(grover, different_targets_return_valid_indices) {
    for (int target = 0; target < 4; target++) {
        int result = Grover(2, target).run();
        cr_assert(result >= 0 && result < 4,
                  "Grover(2,%d) must return an index in [0,4), got %d", target, result);
    }
}

/* ── DeutschJozsa ──────────────────────────────────────────────────────────── */

Test(deutsch_jozsa, constant_zero_oracle_returns_constant) {
    DeutschJozsa dj(1, [](QuantumState&) {});
    cr_assert_eq(dj.run(), FunctionType::CONSTANT);
}

Test(deutsch_jozsa, constant_one_oracle_returns_constant) {
    DeutschJozsa dj(1, [](QuantumState& s) { s.xGate(1); });
    cr_assert_eq(dj.run(), FunctionType::CONSTANT);
}

Test(deutsch_jozsa, balanced_oracle_n1_returns_balanced) {
    DeutschJozsa dj(1, [](QuantumState& s) { s.cnotGate(0, 1); });
    cr_assert_eq(dj.run(), FunctionType::BALANCED);
}

Test(deutsch_jozsa, n2_constant_oracle_returns_constant) {
    DeutschJozsa dj(2, [](QuantumState&) {});
    cr_assert_eq(dj.run(), FunctionType::CONSTANT);
}

Test(deutsch_jozsa, n2_balanced_oracle_returns_balanced) {
    DeutschJozsa dj(2, [](QuantumState& s) { s.cnotGate(0, 2); });
    cr_assert_eq(dj.run(), FunctionType::BALANCED);
}

/* ── swap gate in parser ───────────────────────────────────────────────────── */

Test(quantum_circuit, swap_circuit_swaps_qubits) {
    QuantumCircuit c(2);
    c.load("tests/functional/circuits/swap_circuit.qasm");

    cr_assert_eq(c.measure(), 2);
}

/* ── ccx gate in parser ────────────────────────────────────────────────────── */

Test(quantum_circuit, ccx_circuit_flips_target_when_both_controls_set) {
    QuantumCircuit c(3);
    c.load("tests/functional/circuits/ccx_circuit.qasm");

    cr_assert_eq(c.measure(), 7);
}

/* ── draw ──────────────────────────────────────────────────────────────────── */

Test(quantum_circuit, draw_does_not_crash_on_empty_circuit) {
    QuantumCircuit c(2);
    c.draw();
}

Test(quantum_circuit, draw_outputs_one_line_per_qubit, .init = cr_redirect_stdout) {
    QuantumCircuit c(2);
    c.load("tests/functional/circuits/bell_state.qasm");
    c.draw();
    fflush(stdout);

    FILE* f = cr_get_redirected_stdout();
    rewind(f);
    char buf[1024] = {};
    fread(buf, 1, sizeof(buf) - 1, f);
    std::string out(buf);

    cr_assert(out.find("q[0]") != std::string::npos, "draw() must contain q[0]");
    cr_assert(out.find("q[1]") != std::string::npos, "draw() must contain q[1]");
    cr_assert(out.find("H")    != std::string::npos, "draw() must show H gate");
}

/* ── printState ────────────────────────────────────────────────────────────── */

Test(quantum_circuit, print_state_does_not_crash) {
    QuantumCircuit c(1);
    c.printState();
}

Test(quantum_circuit, print_state_shows_excited_state, .init = cr_redirect_stdout) {
    QuantumCircuit c(1);
    c.load("tests/functional/circuits/x_gate.qasm");
    c.printState();
    fflush(stdout);

    FILE* f = cr_get_redirected_stdout();
    rewind(f);
    char buf[512] = {};
    fread(buf, 1, sizeof(buf) - 1, f);
    std::string out(buf);

    cr_assert(out.find("100") != std::string::npos, "printState must show 100%% probability");
    cr_assert(out.find("1")   != std::string::npos, "printState must show binary '1'");
}

/* ── Simon's algorithm ─────────────────────────────────────────────────────── */

Test(simon, detects_period_s1_n2) {
    Simon simon;
    int s = simon.run(2, [](int x) { return x >> 1; });

    cr_assert_eq(s, 1);
}

Test(simon, detects_period_s2_n2) {
    Simon simon;
    int s = simon.run(2, [](int x) { return x & 1; });

    cr_assert_eq(s, 2);
}

Test(simon, detects_period_s1_n3) {
    Simon simon;
    int s = simon.run(3, [](int x) { return x >> 1; });

    cr_assert_eq(s, 1);
}

Test(simon, detects_period_s2_n3) {
    Simon simon;
    /* secret s=2=010: f(x)=f(x^2), oracle maps each pair to its minimum */
    int s = simon.run(3, [](int x) { int y = x ^ 2; return x < y ? x : y; });

    cr_assert_eq(s, 2);
}

/* ── DensityMatrix: initialisation ────────────────────────────────────────── */

Test(density_matrix, init_ground_state_diagonal) {
    DensityMatrix dm;
    dm.initialize(2);

    cr_assert_float_eq(dm.getDiagonal(0), 1.0, EPS);
    cr_assert_float_eq(dm.getDiagonal(1), 0.0, EPS);
    cr_assert_float_eq(dm.getDiagonal(2), 0.0, EPS);
    cr_assert_float_eq(dm.getDiagonal(3), 0.0, EPS);
}

Test(density_matrix, throws_above_qubit_limit) {
    DensityMatrix dm;
    cr_assert_throw(dm.initialize(DM_MAX_QUBITS + 1), std::runtime_error);
}

/* ── DensityMatrix: pure-state gates ──────────────────────────────────────── */

Test(density_matrix, x_gate_flips_to_excited) {
    DensityMatrix dm;
    dm.initialize(1);
    dm.xGate(0);

    cr_assert_float_eq(dm.getDiagonal(0), 0.0, EPS);
    cr_assert_float_eq(dm.getDiagonal(1), 1.0, EPS);
}

Test(density_matrix, x_gate_twice_is_identity) {
    DensityMatrix dm;
    dm.initialize(1);
    dm.xGate(0);
    dm.xGate(0);

    cr_assert_float_eq(dm.getDiagonal(0), 1.0, EPS);
    cr_assert_float_eq(dm.getDiagonal(1), 0.0, EPS);
}

Test(density_matrix, h_gate_creates_equal_superposition) {
    DensityMatrix dm;
    dm.initialize(1);
    dm.hGate(0);

    cr_assert_float_eq(dm.getDiagonal(0), 0.5, EPS);
    cr_assert_float_eq(dm.getDiagonal(1), 0.5, EPS);
}

Test(density_matrix, h_gate_twice_is_identity) {
    DensityMatrix dm;
    dm.initialize(1);
    dm.hGate(0);
    dm.hGate(0);

    cr_assert_float_eq(dm.getDiagonal(0), 1.0, EPS);
    cr_assert_float_eq(dm.getDiagonal(1), 0.0, EPS);
}

Test(density_matrix, ry_pi_flips_state) {
    /* Ry(π)|0⟩ = |1⟩ — tests gateRight with asymmetric m01≠m10 */
    DensityMatrix dm;
    dm.initialize(1);
    dm.ryGate(M_PI, 0);

    cr_assert_float_eq(dm.getDiagonal(0), 0.0, EPS);
    cr_assert_float_eq(dm.getDiagonal(1), 1.0, EPS);
}

Test(density_matrix, cnot_creates_bell_state) {
    DensityMatrix dm;
    dm.initialize(2);
    dm.hGate(0);
    dm.cnotGate(0, 1);

    cr_assert_float_eq(dm.getDiagonal(0), 0.5, EPS);
    cr_assert_float_eq(dm.getDiagonal(1), 0.0, EPS);
    cr_assert_float_eq(dm.getDiagonal(2), 0.0, EPS);
    cr_assert_float_eq(dm.getDiagonal(3), 0.5, EPS);
}

/* ── DensityMatrix: noise channels ────────────────────────────────────────── */

Test(density_matrix, bit_flip_correct_probs) {
    /* |0⟩ + bit-flip(p=0.20) → 80% |0⟩, 20% |1⟩ */
    DensityMatrix dm;
    dm.initialize(1);
    dm.applyBitFlip(0, 0.20);

    cr_assert_float_eq(dm.getDiagonal(0), 0.80, EPS);
    cr_assert_float_eq(dm.getDiagonal(1), 0.20, EPS);
}

Test(density_matrix, bit_flip_trace_preserved) {
    DensityMatrix dm;
    dm.initialize(1);
    dm.applyBitFlip(0, 0.35);

    cr_assert_float_eq(dm.getDiagonal(0) + dm.getDiagonal(1), 1.0, EPS);
}

Test(density_matrix, phase_flip_full_dephasing) {
    /* |+⟩ + phase-flip(p=1) → |−⟩⟨−|, diagonal [0.5, 0.5] */
    DensityMatrix dm;
    dm.initialize(1);
    dm.hGate(0);
    dm.applyPhaseFlip(0, 1.0);

    cr_assert_float_eq(dm.getDiagonal(0), 0.5, EPS);
    cr_assert_float_eq(dm.getDiagonal(1), 0.5, EPS);
}

Test(density_matrix, phase_flip_trace_preserved) {
    DensityMatrix dm;
    dm.initialize(1);
    dm.hGate(0);
    dm.applyPhaseFlip(0, 0.40);

    cr_assert_float_eq(dm.getDiagonal(0) + dm.getDiagonal(1), 1.0, EPS);
}

Test(density_matrix, depolarizing_correct_probs) {
    /* depolarizing on |0⟩: diag[0] = 1 - 2p/3, diag[1] = 2p/3 */
    DensityMatrix dm;
    dm.initialize(1);
    double p = 0.30;
    dm.applyDepolarizing(0, p);

    cr_assert_float_eq(dm.getDiagonal(0), 1.0 - 2.0 * p / 3.0, EPS);
    cr_assert_float_eq(dm.getDiagonal(1), 2.0 * p / 3.0,       EPS);
}

Test(density_matrix, depolarizing_trace_preserved) {
    DensityMatrix dm;
    dm.initialize(1);
    dm.applyDepolarizing(0, 0.45);

    cr_assert_float_eq(dm.getDiagonal(0) + dm.getDiagonal(1), 1.0, EPS);
}

Test(density_matrix, depolarizing_multiple_gates_trace_preserved) {
    /* Bell state + two depolarizing channels */
    DensityMatrix dm;
    dm.initialize(2);
    dm.hGate(0);
    dm.cnotGate(0, 1);
    dm.applyDepolarizing(0, 0.10);
    dm.applyDepolarizing(1, 0.10);

    double trace = dm.getDiagonal(0) + dm.getDiagonal(1)
                 + dm.getDiagonal(2) + dm.getDiagonal(3);
    cr_assert_float_eq(trace, 1.0, EPS);
}

/* ── DensityMatrix: QuantumCircuit integration ─────────────────────────────── */

Test(density_matrix, circuit_valid_noise_throws_above_limit) {
    QuantumCircuit circuit(13);
    cr_assert_throw(circuit.validNoise(0.01), std::runtime_error);
}

Test(density_matrix, circuit_noise_measure_in_range) {
    QuantumCircuit circuit(2);
    circuit.validNoise(0.01);
    circuit.load("tests/functional/circuits/bell_state.qasm");
    int result = circuit.measure();
    cr_assert_geq(result, 0);
    cr_assert_lt(result, 4);
}

Test(simon, detects_period_s4_n3) {
    Simon simon;
    /* secret s=4=100: f(x)=f(x^4) */
    int s = simon.run(3, [](int x) { int y = x ^ 4; return x < y ? x : y; });

    cr_assert_eq(s, 4);
}

Test(simon, detects_period_s1_n4) {
    Simon simon;
    /* secret s=1: f(x)=x>>1 works for any n since the LSB is the period */
    int s = simon.run(4, [](int x) { return x >> 1; });

    cr_assert_eq(s, 1);
}

/* ── Shor's algorithm ──────────────────────────────────────────────────────── */

Test(shor, modular_oracle_does_not_crash) {
    Shor shor;
    cr_assert_no_throw(shor.run(15), std::exception);
}

Test(shor, run_returns_valid_factors_of_15) {
    Shor shor;
    auto [p, q] = shor.run(15);
    if (p != -1) {
        cr_assert_eq(p * q, 15, "p * q doit être 15");
        cr_assert(p > 1 && q > 1, "les deux facteurs doivent être non-triviaux");
    }
}

Test(shor, run_returns_valid_factors_of_21) {
    Shor shor;
    auto [p, q] = shor.run(21);
    if (p != -1) {
        cr_assert_eq(p * q, 21, "p * q doit être 21");
        cr_assert(p > 1 && q > 1, "les deux facteurs doivent être non-triviaux");
    }
}

/* ── VQE ───────────────────────────────────────────────────────────────────── */

Test(vqe, constructor_does_not_crash) {
    VQE vqe(1, {{1.0, {PauliOp::Z}}});
}

Test(vqe, get_optimal_params_is_empty_before_run) {
    VQE vqe(1, {{1.0, {PauliOp::Z}}});

    cr_assert_eq(vqe.getOptimalParams().size(), 0);
}

Test(vqe, run_converges_toward_ground_state_single_z) {
    VQE vqe(1, {{1.0, {PauliOp::Z}}});

    double energy = vqe.run();
    cr_assert_lt(energy, -0.9, "VQE doit converger vers E ≈ −1 pour H = Z₀");
}

Test(vqe, get_optimal_params_has_correct_size_after_run) {
    VQE vqe(2, {{1.0, {PauliOp::Z, PauliOp::Z}}});
    vqe.run();

    cr_assert_eq(vqe.getOptimalParams().size(), 2u);
}

Test(vqe, run_converges_toward_ground_state_two_qubit_zz) {
    VQE vqe(2, {{1.0, {PauliOp::Z, PauliOp::Z}}});

    double energy = vqe.run();
    cr_assert_lt(energy, -0.9, "VQE doit converger vers E ≈ −1 pour H = Z₀⊗Z₁");
}

Test(vqe, run_with_identity_term_returns_coefficient) {
    VQE vqe(1, {{2.5, {PauliOp::I}}});

    double energy = vqe.run();
    cr_assert_float_eq(energy, 2.5, 1e-9);
}

/* ── Shor — failure path ────────────────────────────────────────────────────── */

Test(shor, run_returns_failure_for_prime_input) {
    Shor shor;
    auto [p, q] = shor.run(5);
    cr_assert_eq(p, -1);
    cr_assert_eq(q, -1);
}

/* ── draw — swap and ccx symbols ────────────────────────────────────────────── */

Test(quantum_circuit, draw_includes_swap_symbol, .init = cr_redirect_stdout) {
    QuantumCircuit c(2);
    c.load("tests/functional/circuits/swap_circuit.qasm");
    c.draw();
    fflush(stdout);

    FILE* f = cr_get_redirected_stdout();
    rewind(f);
    char buf[1024] = {};
    fread(buf, 1, sizeof(buf) - 1, f);
    std::string out(buf);

    cr_assert(out.find("×") != std::string::npos, "draw() must show × for swap gate");
}

Test(quantum_circuit, draw_includes_ccx_symbols, .init = cr_redirect_stdout) {
    QuantumCircuit c(3);
    c.load("tests/functional/circuits/ccx_circuit.qasm");
    c.draw();
    fflush(stdout);

    FILE* f = cr_get_redirected_stdout();
    rewind(f);
    char buf[1024] = {};
    fread(buf, 1, sizeof(buf) - 1, f);
    std::string out(buf);

    cr_assert(out.find("●") != std::string::npos, "draw() must show ● for ccx control qubits");
}

/* ── draw — non-adjacent cx produces │ pipe between qubits ─────────────────── */

Test(quantum_circuit, draw_pipe_for_nonadjacent_cx, .init = cr_redirect_stdout) {
    QuantumCircuit c(3);
    c.load("tests/functional/circuits/cx_nonadjacent.qasm");
    c.draw();
    fflush(stdout);

    FILE* f = cr_get_redirected_stdout();
    rewind(f);
    char buf[1024] = {};
    fread(buf, 1, sizeof(buf) - 1, f);
    std::string out(buf);

    cr_assert(out.find("│") != std::string::npos, "draw() must show │ for qubits between cx operands");
    cr_assert(out.find("●") != std::string::npos, "draw() must show ● for cx control qubit");
}
