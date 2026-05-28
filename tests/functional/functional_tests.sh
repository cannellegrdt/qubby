#!/usr/bin/env bash
#
# Project: qubby
# File name: functional_tests.sh
# Author: Cannelle Gourdet - lankley
# File description: Functional tests for the qubby CLI.
#                   Each test runs the binary with a .qasm circuit file and
#                   checks exit code and/or output.

set -euo pipefail

BINARY="./qubby"
CIRCUITS="tests/functional/circuits"
PASS=0
FAIL=0

RED='\033[0;31m'
GREEN='\033[0;32m'
RESET='\033[0m'

pass() { echo -e "  ${GREEN}[PASS]${RESET} $1"; PASS=$((PASS + 1)); }
fail() { echo -e "  ${RED}[FAIL]${RESET} $1"; FAIL=$((FAIL + 1)); }

# run_test <name> <expected_exit> <pattern> <cmd...>
#   pattern: fixed string that stdout+stderr must contain, or "" to skip.
#            Pipe-separate alternatives for OR matching: "foo|bar"
run_test() {
    local name="$1"
    local expected_exit="$2"
    local pattern="$3"
    shift 3

    local output
    output=$("$@" 2>&1) && actual_exit=0 || actual_exit=$?

    if [ "$actual_exit" -ne "$expected_exit" ]; then
        fail "$name — exit $actual_exit (expected $expected_exit)"
        return
    fi

    if [ -n "$pattern" ] && ! echo "$output" | grep -qE "$pattern"; then
        fail "$name — output «$output» does not match «$pattern»"
        return
    fi

    pass "$name"
}


echo "── Functional tests ─────────────────────────────────────────────────────"

# 1. X gate: |0⟩ → |1⟩  →  measurement is always 1
run_test "x_gate_deterministic_measure" \
    0 "Measurement: 1" \
    "$BINARY" "$CIRCUITS/x_gate.qasm" 1

# 2. H²=I: H then H brings state back to |0⟩  →  measurement is always 0
run_test "hh_identity_measure" \
    0 "Measurement: 0" \
    "$BINARY" "$CIRCUITS/hh_identity.qasm" 1

# 3. Comments: lines starting with // must be ignored, only X gate runs
run_test "comments_ignored" \
    0 "Measurement: 1" \
    "$BINARY" "$CIRCUITS/comments.qasm" 1

# 4. Z|0⟩ = |0⟩ : Z on ground state does not change measurement
run_test "z_on_ground_state_measures_zero" \
    0 "Measurement: 0" \
    "$BINARY" "$CIRCUITS/z_gate.qasm" 1

# 5. X then Z: phase flip on |1⟩ = -|1⟩, measurement probability unchanged
run_test "xz_excited_state_still_measures_one" \
    0 "Measurement: 1" \
    "$BINARY" "$CIRCUITS/xz_gate.qasm" 1

# 6. Bell state: must run without error and return 0 or 3
run_test "bell_state_valid_output" \
    0 "Measurement: 0|Measurement: 3" \
    "$BINARY" "$CIRCUITS/bell_state.qasm" 2

# 7. Missing file must exit non-zero
run_test "missing_file_exits_nonzero" \
    1 "Error:" \
    "$BINARY" "$CIRCUITS/nonexistent.qasm" 1

# 8. Unknown gate must exit non-zero
run_test "unknown_gate_exits_nonzero" \
    1 "Error:" \
    "$BINARY" "$CIRCUITS/invalid_gate.qasm" 1

# 9. Wrong number of arguments
run_test "no_args_exits_nonzero" \
    1 "" \
    "$BINARY"

# 10. Rx(π)|0⟩ = -i|1⟩ → probability 1, measurement is always 1
run_test "rx_pi_measures_one" \
    0 "Measurement: 1" \
    "$BINARY" "$CIRCUITS/rx_pi.qasm" 1

# 11. Ry(π)|0⟩ = |1⟩ → measurement is always 1
run_test "ry_pi_measures_one" \
    0 "Measurement: 1" \
    "$BINARY" "$CIRCUITS/ry_pi.qasm" 1

# 12. Rz(π/2)|0⟩: phase change only, measurement is always 0
run_test "rz_on_ground_state_measures_zero" \
    0 "Measurement: 0" \
    "$BINARY" "$CIRCUITS/rz_ground.qasm" 1

# 13. Block comment /* */ is ignored: only the X gate runs → measurement = 1
run_test "block_comment_ignored" \
    0 "Measurement: 1" \
    "$BINARY" "$CIRCUITS/block_comment.qasm" 1

# 14. qreg reinitialises the state, then X runs → measurement = 1
run_test "qreg_reinitializes_state" \
    0 "Measurement: 1" \
    "$BINARY" "$CIRCUITS/qreg.qasm" 1

# 15. Symbolic angle rx(pi) → Rx(π)|0⟩ = -i|1⟩ → measurement = 1
run_test "rx_symbolic_pi_measures_one" \
    0 "Measurement: 1" \
    "$BINARY" "$CIRCUITS/rx_symbolic_pi.qasm" 1

# 16. Symbolic angle ry(pi/2) → equal superposition, result is 0 or 1
run_test "ry_symbolic_pi_div2_valid_output" \
    0 "Measurement: 0|Measurement: 1" \
    "$BINARY" "$CIRCUITS/ry_symbolic_pi_div2.qasm" 1

# 17. Multi-line block comment /* ... */ is ignored, X gate runs → measurement = 1
run_test "multiline_block_comment_ignored" \
    0 "Measurement: 1" \
    "$BINARY" "$CIRCUITS/multiline_block_comment.qasm" 1

# 18. Symbolic angle rx(pi*2) → Rx(2π)|0⟩ = -|0⟩ → measurement = 0
run_test "rx_symbolic_pi_mult_measures_zero" \
    0 "Measurement: 0" \
    "$BINARY" "$CIRCUITS/rx_symbolic_pi_mult.qasm" 1

# 19. Malformed rx angle (missing parentheses) → parsing error
run_test "malformed_rx_angle_exits_nonzero" \
    1 "Error:" \
    "$BINARY" "$CIRCUITS/malformed_rx_angle.qasm" 1

# 20. --noise produces a valid measurement (0 or 1 for a 1-qubit circuit)
run_test "noise_x_gate_valid_measurement" \
    0 "Measurement: 0|Measurement: 1" \
    "$BINARY" "$CIRCUITS/x_gate.qasm" 1 --noise

# 21. --noise on a Bell state (2 qubits) → valid measurement in 0..3
run_test "noise_bell_state_valid_measurement" \
    0 "Measurement: 0|Measurement: 1|Measurement: 2|Measurement: 3" \
    "$BINARY" "$CIRCUITS/bell_state.qasm" 2 --noise

# 22. --noise-rate=0.50 is accepted (high rate, but valid)
run_test "noise_rate_high_accepted" \
    0 "Measurement:" \
    "$BINARY" "$CIRCUITS/x_gate.qasm" 1 --noise --noise-rate=0.50

# 23. --noise --print displays percentages (output contains %)
run_test "noise_print_shows_percentages" \
    0 "%" \
    "$BINARY" "$CIRCUITS/bell_state.qasm" 2 --noise --print

# 24. --noise with >12 qubits → exit 1 + error message
run_test "noise_too_many_qubits_exits_nonzero" \
    1 "Error:" \
    "$BINARY" "$CIRCUITS/x_gate.qasm" 13 --noise

# 25. --noise --ascii draws the circuit wires (q[0]: ...) and produces a measurement
run_test "noise_ascii_draw_shows_wires" \
    0 "q\[0\]:" \
    "$BINARY" "$CIRCUITS/bell_state.qasm" 2 --noise --ascii

run_test "noise_ascii_draw_ends_with_measurement" \
    0 "Measurement:" \
    "$BINARY" "$CIRCUITS/bell_state.qasm" 2 --noise --ascii

# 27. --bloch=0 on x_gate (|1⟩): output contains |r| and south pole info
run_test "bloch_x_gate_shows_r_equal_one" \
    0 "\|r\| = 1" \
    "$BINARY" "$CIRCUITS/x_gate.qasm" 1 --bloch=0

# 28. --bloch=0 on ground state: Bloch z component is positive (north pole)
run_test "bloch_ground_state_shows_output" \
    0 "Bloch sphere" \
    "$BINARY" "$CIRCUITS/hh_identity.qasm" 1 --bloch=0

# 29. --bloch with --noise active: prints a warning, still exits 0
run_test "bloch_with_noise_prints_warning" \
    0 "Warning:" \
    "$BINARY" "$CIRCUITS/x_gate.qasm" 1 --noise --bloch=0

# 30. --bloch-svg + --bloch-file creates an SVG file containing '<svg'
run_test "bloch_svg_creates_valid_file" \
    0 "Saved Bloch sphere SVG" \
    "$BINARY" "$CIRCUITS/x_gate.qasm" 1 --bloch-svg=0 --bloch-file=/tmp/qubby_test.svg
if grep -q "<svg" /tmp/qubby_test.svg 2>/dev/null; then
    pass "bloch_svg_file_contains_svg_element"
else
    fail "bloch_svg_file_contains_svg_element"
fi

# 32. --export-qiskit produces a Python file with QuantumCircuit
run_test "export_qiskit_contains_circuit" \
    0 "Exported Qiskit" \
    "$BINARY" "$CIRCUITS/bell_state.qasm" 2 --export-qiskit=/tmp/qubby_bell.py
if grep -q "QuantumCircuit" /tmp/qubby_bell.py 2>/dev/null; then
    pass "export_qiskit_file_has_QuantumCircuit"
else
    fail "export_qiskit_file_has_QuantumCircuit"
fi

# 34. --export-cirq produces a Python file with 'import cirq'
run_test "export_cirq_contains_circuit" \
    0 "Exported Cirq" \
    "$BINARY" "$CIRCUITS/bell_state.qasm" 2 --export-cirq=/tmp/qubby_bell_cirq.py
if grep -q "import cirq" /tmp/qubby_bell_cirq.py 2>/dev/null; then
    pass "export_cirq_file_has_import_cirq"
else
    fail "export_cirq_file_has_import_cirq"
fi

# 36. --repl with EOF (empty stdin) exits cleanly (code 0)
run_test "repl_exits_cleanly_on_eof" \
    0 "Bye" \
    bash -c "echo '' | $BINARY --repl"

# 37. --repl accepts qreg + h + measure and prints state
run_test "repl_executes_instructions" \
    0 "%" \
    bash -c "printf 'h q[0]\n.quit\n' | $BINARY --repl 1"

# 38. --repl .help prints command list
run_test "repl_help_shows_commands" \
    0 "\.quit" \
    bash -c "printf '.help\n.quit\n' | $BINARY --repl"

# 39. --repl .measure after X gate prints 1
run_test "repl_measure_after_x_returns_one" \
    0 "Measured: 1" \
    bash -c "printf 'x q[0]\n.measure\n.quit\n' | $BINARY --repl 1"

# 40. --repl .reset clears state back to ground
run_test "repl_reset_clears_state" \
    0 "Measured: 0" \
    bash -c "printf 'x q[0]\n.reset\n.measure\n.quit\n' | $BINARY --repl 1"

# 41. --repl .ascii draws circuit wires
run_test "repl_ascii_draws_wires" \
    0 "q\[0\]:" \
    bash -c "printf 'h q[0]\n.ascii\n.quit\n' | $BINARY --repl 1"

# 42. --repl .export-qiskit creates a file
run_test "repl_export_qiskit_creates_file" \
    0 "Exported Qiskit" \
    bash -c "printf 'h q[0]\ncx q[0] q[1]\n.export-qiskit /tmp/repl_export.py\n.quit\n' | $BINARY --repl 2"

echo ""
echo "  Results: ${PASS} passed, ${FAIL} failed"
echo ""
[ "$FAIL" -eq 0 ]
