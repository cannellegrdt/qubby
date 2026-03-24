/*
 * Project: qubby
 * File name: test_sycl.cpp
 * Author: Cannelle Gourdet - lankley
 * File description: Criterion unit tests for QuantumState::applyGateSYCL.
 *
 * All tests are guarded by #ifdef SYCL_ENABLED — this file compiles cleanly
 * in a standard g++ build (no tests emitted) and is fully exercised when
 * built with: make sycl_tests  (uses icpx -DSYCL_ENABLED)
 *
 * Strategy: for each gate, verify both the raw amplitude values AND that the
 * SYCL kernel produces bit-for-bit identical results to the CPU applyGate path.
 */

#ifdef SYCL_ENABLED

#include <criterion/criterion.h>
#include "QuantumState.hpp"
#include <cmath>

static constexpr double EPS = 1e-9;

static void assert_sycl_matches_cpu(int numQubits, int targetQubit, Matrix gate) {
    QuantumState cpu, sycl_s;
    cpu.initialize(numQubits);
    sycl_s.initialize(numQubits);

    cpu.applyGate(targetQubit, gate);
    sycl_s.applyGateSYCL(targetQubit, gate);

    int size = 1 << numQubits;
    for (int i = 0; i < size; i++) {
        cr_assert_float_eq(sycl_s.getAmplitude(i).real(),
                           cpu.getAmplitude(i).real(), EPS,
                           "amplitude[%d].real mismatch", i);
        cr_assert_float_eq(sycl_s.getAmplitude(i).imag(),
                           cpu.getAmplitude(i).imag(), EPS,
                           "amplitude[%d].imag mismatch", i);
    }
}

/* ── X gate ─────────────────────────────────────────────────────────────────── */

Test(sycl_x_gate, flips_ground_to_excited) {
    QuantumState s;
    s.initialize(1);
    Matrix x = {0, 1, 1, 0};
    s.applyGateSYCL(0, x);

    cr_assert_float_eq(s.getAmplitude(0).real(), 0.0, EPS);
    cr_assert_float_eq(s.getAmplitude(1).real(), 1.0, EPS);
}

Test(sycl_x_gate, applied_twice_is_identity) {
    QuantumState s;
    s.initialize(1);
    Matrix x = {0, 1, 1, 0};
    s.applyGateSYCL(0, x);
    s.applyGateSYCL(0, x);

    cr_assert_float_eq(s.getAmplitude(0).real(), 1.0, EPS);
    cr_assert_float_eq(s.getAmplitude(1).real(), 0.0, EPS);
}

Test(sycl_x_gate, targets_correct_qubit_in_two_qubit_system) {
    QuantumState s;
    s.initialize(2);
    Matrix x = {0, 1, 1, 0};
    s.applyGateSYCL(1, x);

    cr_assert_float_eq(s.getAmplitude(0).real(), 0.0, EPS);
    cr_assert_float_eq(s.getAmplitude(2).real(), 1.0, EPS);
}

Test(sycl_x_gate, matches_cpu_path_single_qubit) {
    Matrix x = {0, 1, 1, 0};
    assert_sycl_matches_cpu(1, 0, x);
}

Test(sycl_x_gate, matches_cpu_path_three_qubits_middle_qubit) {
    Matrix x = {0, 1, 1, 0};
    assert_sycl_matches_cpu(3, 1, x);
}

/* ── H gate ─────────────────────────────────────────────────────────────────── */

Test(sycl_h_gate, creates_equal_superposition) {
    QuantumState s;
    s.initialize(1);
    double v = 1.0 / std::sqrt(2.0);
    Matrix h = {v, v, v, -v};
    s.applyGateSYCL(0, h);

    cr_assert_float_eq(s.getAmplitude(0).real(), v, EPS);
    cr_assert_float_eq(s.getAmplitude(1).real(), v, EPS);
}

Test(sycl_h_gate, applied_twice_is_identity) {
    QuantumState s;
    s.initialize(1);
    double v = 1.0 / std::sqrt(2.0);
    Matrix h = {v, v, v, -v};
    s.applyGateSYCL(0, h);
    s.applyGateSYCL(0, h);

    cr_assert_float_eq(s.getAmplitude(0).real(), 1.0, EPS);
    cr_assert_float_eq(s.getAmplitude(1).real(), 0.0, EPS);
}

Test(sycl_h_gate, matches_cpu_path_single_qubit) {
    double v = 1.0 / std::sqrt(2.0);
    Matrix h = {v, v, v, -v};
    assert_sycl_matches_cpu(1, 0, h);
}

Test(sycl_h_gate, matches_cpu_path_two_qubits_second_qubit) {
    double v = 1.0 / std::sqrt(2.0);
    Matrix h = {v, v, v, -v};
    assert_sycl_matches_cpu(2, 1, h);
}

/* ── Z gate ─────────────────────────────────────────────────────────────────── */

Test(sycl_z_gate, leaves_ground_state_unchanged) {
    QuantumState s;
    s.initialize(1);
    Matrix z = {1, 0, 0, -1};
    s.applyGateSYCL(0, z);

    cr_assert_float_eq(s.getAmplitude(0).real(),  1.0, EPS);
    cr_assert_float_eq(s.getAmplitude(1).real(),  0.0, EPS);
}

Test(sycl_z_gate, negates_excited_state_amplitude) {
    QuantumState s;
    s.initialize(1);
    Matrix x = {0, 1, 1, 0};
    Matrix z = {1, 0, 0, -1};
    s.applyGateSYCL(0, x);
    s.applyGateSYCL(0, z);

    cr_assert_float_eq(s.getAmplitude(0).real(),  0.0, EPS);
    cr_assert_float_eq(s.getAmplitude(1).real(), -1.0, EPS);
}

Test(sycl_z_gate, hzh_equals_x) {
    QuantumState s;
    s.initialize(1);
    double v = 1.0 / std::sqrt(2.0);
    Matrix h = {v, v, v, -v};
    Matrix z = {1, 0, 0, -1};
    s.applyGateSYCL(0, h);
    s.applyGateSYCL(0, z);
    s.applyGateSYCL(0, h);

    cr_assert_float_eq(s.getAmplitude(0).real(), 0.0, EPS);
    cr_assert_float_eq(s.getAmplitude(1).real(), 1.0, EPS);
}

Test(sycl_z_gate, matches_cpu_path_single_qubit) {
    Matrix z = {1, 0, 0, -1};
    assert_sycl_matches_cpu(1, 0, z);
}

/* ── cross-path: Bell state with SYCL H + CPU CNOT ─────────────────────────── */

Test(sycl_bell_state, h_via_sycl_cnot_via_cpu) {
    QuantumState s;
    s.initialize(2);
    double v = 1.0 / std::sqrt(2.0);
    Matrix h = {v, v, v, -v};
    s.applyGateSYCL(0, h);
    s.cnotGate(0, 1);

    cr_assert_float_eq(s.getAmplitude(0).real(), v, EPS);
    cr_assert_float_eq(s.getAmplitude(3).real(), v, EPS);
    cr_assert_float_eq(s.getAmplitude(1).real(), 0.0, EPS);
    cr_assert_float_eq(s.getAmplitude(2).real(), 0.0, EPS);
}

#endif /* SYCL_ENABLED */
