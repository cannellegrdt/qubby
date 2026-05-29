# qubby

A quantum circuit simulator written in C++20, capable of simulating up to 25 qubits. It runs OpenQASM 2.0 circuit files from the command line and includes hardware acceleration via OpenMP (CPU multithreading) and SYCL via Intel oneAPI.

## Features

- **State vector simulation** — contiguous `std::complex<double>` vector of size 2ⁿ, up to 25 qubits (512 MB)
- **Bitwise gate application** — in-place 2×2 matrix transforms using bit-index masking; no full 2ⁿ×2ⁿ matrices
- **OpenQASM 2.0 parser** — load and run `.qasm` circuit files with register declarations, parametric gates, and block comments
- **Hardware acceleration** — OpenMP for CPU parallelism; SYCL kernel for Intel GPUs via oneAPI (optional)
- **Rich gate set** — X, H, Z, Rx/Ry/Rz, CNOT, SWAP, Toffoli, controlled-Rz
- **Quantum algorithms** — Deutsch-Jozsa, Grover, QFT, Shor, Simon, VQE
- **Noise simulation** — depolarizing, bit-flip and phase-flip channels on a density matrix backend
- **Bloch sphere** — ASCII terminal view + SVG export of any single-qubit state
- **Cross-platform export** — generate equivalent Qiskit or Cirq Python scripts from a loaded circuit
- **Interactive REPL** — type OpenQASM instructions in real time and watch the state evolve live

## Gate set

| Gate | Description |
|------|-------------|
| `x` | Pauli-X (NOT) |
| `h` | Hadamard (superposition) |
| `z` | Pauli-Z (phase flip) |
| `rx(θ)` | Rotation around X-axis |
| `ry(θ)` | Rotation around Y-axis |
| `rz(θ)` | Rotation around Z-axis |
| `cx` | Controlled-NOT (CNOT) |
| `swap` | SWAP |
| `ccx` | Toffoli (CCX) |

## Quantum algorithms

| Algorithm | Complexity | Description |
|-----------|------------|-------------|
| Deutsch-Jozsa | O(1) queries | Determines if a function is constant or balanced |
| Grover | O(√N) | Amplitude amplification for unstructured search |
| QFT | O(n²) | Quantum Fourier Transform — core subroutine for Shor |
| Shor | O((log N)³) | Integer factorisation in polynomial time |
| Simon | O(n) queries | Finds a hidden period exponentially faster than classical |
| VQE | Hybrid | Variational Quantum Eigensolver for ground-state energy estimation |

## Requirements

- C++20 compiler (`g++`)
- OpenMP
- [Criterion](https://github.com/Snaipe/Criterion) (unit tests only)
- `lcov` + `genhtml` (coverage only)
- Intel oneAPI Base Toolkit with `icpx` (SYCL build only)

## Build

```bash
# Standard build (OpenMP)
make

# Clean rebuild
make re

# Build with Intel oneAPI SYCL (requires icpx)
make sycl
```

## Usage

```
./qubby <circuit.qasm> <num_qubits> [OPTIONS]
./qubby --repl [num_qubits]

Arguments:
  <circuit.qasm>           Path to an OpenQASM 2.0 circuit file
  <num_qubits>             Number of qubits in the circuit

Options:
  --repl                   Start the interactive OpenQASM REPL
  --noise                  Enable noise simulation (density matrix, ≤12 qubits)
  --noise-rate=<p>         Per-gate depolarizing error rate in [0,1] (default: 0.01)
  --ascii                  Draw an ASCII diagram of the circuit
  --print                  Print the quantum state before measurement
  --bloch=<qubit>          Show Bloch sphere (ASCII) for the given qubit index
  --bloch-svg=<qubit>      Same, and save an SVG (requires --bloch-file)
  --bloch-file=<path>      Output path for the Bloch sphere SVG
  --export-qiskit=<file>   Export circuit as a Qiskit Python script
  --export-cirq=<file>     Export circuit as a Cirq Python script
  -h, --help               Show this help message and exit
```

## Examples

### Bell state — ASCII diagram + probability distribution

```bash
./qubby tests/functional/circuits/bell_state.qasm 2 --ascii --print
```

```
q[0]: ----H--●----
q[1]: -------X----
| 0⟩  00  0.707+0.000i  (50.000%)
| 3⟩  11  0.707+0.000i  (50.000%)
Measurement: 3
```

### Bloch sphere — visualise a single-qubit state

```bash
# |1⟩ state (X gate) → south pole
./qubby tests/functional/circuits/x_gate.qasm 1 --bloch=0
```

```
  Bloch sphere — q[0]  (xz projection)

       |0⟩
                    |
               ooooo|ooooo
            oooo    |    oooo
          oo        |        oo
        oo          |          oo
       oo           |           oo
      oo            |            oo
      o             |             o
      o             |             o
-x -----------------+----------------- +x
      o             |             o
      o             |             o
      oo            |            oo
       oo           |           oo
        oo          |          oo
          oo        |        oo
            oooo    |    oooo
               ooooo*ooooo
                    |
       |1⟩

  r = (0.0000, -0.0000, -1.0000)   |r| = 1.0000 (pure)
```

```bash
# Save as SVG (open in any browser)
./qubby tests/functional/circuits/x_gate.qasm 1 --bloch-svg=0 --bloch-file=/tmp/bloch.svg
xdg-open /tmp/bloch.svg
```

### Noise simulation — depolarizing channel on a Bell state

```bash
./qubby tests/functional/circuits/bell_state.qasm 2 --noise --noise-rate=0.05 --print
```

```
| 0⟩  00  0.474+0.000i  (47.368%)
| 1⟩  01  0.026+0.000i   (2.632%)
| 2⟩  10  0.026+0.000i   (2.632%)
| 3⟩  11  0.474+0.000i  (47.368%)
Measurement: 0
```

### Export to Qiskit / Cirq

```bash
./qubby tests/functional/circuits/bell_state.qasm 2 --export-qiskit=/tmp/bell.py
cat /tmp/bell.py
```

```python
from qiskit import QuantumCircuit

qc = QuantumCircuit(2)
qc.h(0)
qc.cx(0, 1)

qc.measure_all()
print(qc.draw())
```

```bash
./qubby tests/functional/circuits/bell_state.qasm 2 --export-cirq=/tmp/bell_cirq.py
python3 /tmp/bell_cirq.py   # requires cirq installed
```

### Interactive REPL

```bash
./qubby --repl 2
```

```
qubby REPL  —  type OpenQASM instructions, .help for help, .quit to exit
Initial state: 2 qubit(s) in |0...0⟩

qubby> h q[0]
qubby> cx q[0] q[1]
qubby> .state
| 0⟩  00  0.707+0.000i  (50.000%)
| 3⟩  11  0.707+0.000i  (50.000%)
qubby> .ascii
q[0]: ----H--●----
q[1]: -------X----
qubby> .bloch 0
  Bloch sphere — q[0]  (xz projection)
  r = (0.0000, 0.0000, 0.0000)   |r| = 0.0000 (fully mixed)
qubby> .measure
  Measured: 3 (|11⟩)
qubby> .reset
  Reset to |0...0⟩ (2 qubits)
qubby> .quit
Bye.
```

#### REPL dot-commands

| Command | Description |
|---------|-------------|
| `.state` | Print current amplitude distribution |
| `.measure` | Collapse the state and print the outcome |
| `.reset` | Reinitialise to \|0…0⟩ |
| `.ascii` | Draw the circuit logged so far |
| `.bloch <q>` | Bloch sphere (ASCII) for qubit `q` |
| `.bloch-svg <q> <file>` | Save Bloch sphere SVG for qubit `q` |
| `.export-qiskit <file>` | Export circuit as a Qiskit script |
| `.export-cirq <file>` | Export circuit as a Cirq script |
| `.help` | List all commands |
| `.quit` | Exit the REPL |

### Example OpenQASM circuit

```qasm
qreg q[2];
h q[0];
cx q[0], q[1];
```

## Testing

```bash
# Build and run unit + functional tests
make run_tests

# Unit tests only (153 tests via Criterion)
make unit_tests && ./tests/unit_tests

# Functional tests only (42 shell tests)
make functional_tests

# Generate HTML coverage report
make coverage
# Report available at coverage/html/index.html
```

## Benchmarking

Measures wall-clock time for `applyGate` (OpenMP) vs `applyGateSYCL` (SYCL) across 20–25 qubit counts.

```bash
make benchmark
make benchmark_sycl   # requires icpx
```

## Makefile targets

| Target | Description |
|--------|-------------|
| `all` | Build `qubby` |
| `re` | Clean rebuild |
| `clean` | Remove objects and coverage data |
| `fclean` | `clean` + remove all binaries |
| `unit_tests` | Build unit test binary |
| `functional_tests` | Run functional test suite |
| `run_tests` | Build and run both test suites |
| `coverage` | Build with gcov, run tests, generate HTML report |
| `debug` | Build with ASan/UBSan, no optimisations |
| `sycl` | Build with Intel oneAPI (`icpx -fsycl`) |
| `sycl_tests` | Build and run SYCL unit tests |
| `benchmark` | Build and run OpenMP benchmark |
| `benchmark_sycl` | Build and run OpenMP + SYCL benchmark |
| `format` | Run `clang-format` on all sources |
| `lint` | Run `clang-tidy` on all sources |

## Project structure

```
qubby/
├── src/
│   ├── main.cpp                           # CLI entry point
│   ├── QuantumState.cpp                   # State vector + gate logic
│   ├── QuantumState.sycl.cpp              # SYCL GPU kernel
│   ├── QuantumCircuit.cpp                 # OpenQASM parser and circuit runner
│   ├── BlochSphere.cpp                    # Bloch sphere ASCII + SVG renderer
│   ├── REPL.cpp                           # Interactive OpenQASM REPL
│   ├── DensityMatrix.cpp                  # Density matrix + noise channels
│   ├── DeutschJozsa.cpp                   # Deutsch-Jozsa algorithm (constant vs balanced oracle)
│   ├── Grover.cpp                         # Grover's search via amplitude amplification
│   ├── QuantumFourierTransform.cpp        # QFT — core subroutine for Shor
│   ├── Shor.cpp                           # Shor's integer factorisation algorithm
│   ├── Simon.cpp                          # Simon's hidden-period algorithm
│   └── VariationalQuantumEigensolver.cpp  # VQE hybrid classical/quantum optimiser
├── include/                               # Header files
├── tests/
│   ├── unit/                              # Criterion unit tests (153 tests)
│   └── functional/                        # Shell tests + .qasm circuits (42 tests)
├── benchmarks/
└── Makefile
```

## License

MIT
