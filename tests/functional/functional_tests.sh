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


echo ""
echo "  Results: ${PASS} passed, ${FAIL} failed"
echo ""
[ "$FAIL" -eq 0 ]
