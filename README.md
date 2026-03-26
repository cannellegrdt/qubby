# qubby

A quantum circuit simulator written in C++20, capable of simulating up to 25 qubits. It runs OpenQASM 2.0 circuit files from the command line and includes hardware acceleration via OpenMP (CPU multithreading) and SYCL via Intel oneAPI.

## Features

- **State vector simulation** - contiguous `std::complex<double>` vector of size 2ⁿ, up to 25 qubits (512 MB)
- **Bitwise gate application** - in-place 2×2 matrix transforms using bit-index masking; no full 2ⁿ×2ⁿ matrices
- **OpenQASM 2.0 parser** - load and run `.qasm` circuit files with register declarations, parametric gates, and block comments
- **Hardware acceleration** - OpenMP for CPU parallelism; SYCL kernel for Intel GPUs via oneAPI (optional) - targets any device selected by `sycl::default_selector_v`
- **Rich gate set** - X, H, Z, Rx/Ry/Rz, CNOT, SWAP, Toffoli, controlled-Rz
- **Quantum algorithms** - Deutsch-Jozsa, Grover, QFT, Shor, Simon, VQE
- **Observability** - probability distribution display (`--print`), ASCII circuit drawer (`--ascii`), OpenMP vs SYCL benchmark

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
| QFT | O(n²) | Quantum Fourier Transform - core subroutine for Shor |
| Shor | O((log N)³) | Integer factorisation in polynomial time |
| Simon | O(n) queries | Finds a hidden period exponentially faster than classical |
| VQE | Hybrid | Variational Quantum Eigensolver for ground-state energy estimation |

## Requirements

- C++20 compiler (`g++`)
- OpenMP
- [Criterion](https://github.com/Snaipe/Criterion) (unit tests only)
- `lcov` + `genhtml` (coverage only)
- Intel oneAPI Base Toolkit with `icpx` (SYCL build only - targets any SYCL-compatible device: Intel CPU, Intel GPU, or others via additional plugins)

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

Arguments:
  <circuit.qasm>   Path to an OpenQASM 2.0 circuit file
  <num_qubits>     Number of qubits in the circuit

Options:
  --ascii          Draw an ASCII representation of the circuit
  --print          Print the quantum state before measurement
  -h, --help       Show this help message and exit
```

### Examples

```bash
# Run a Bell state circuit and print the result
./qubby tests/functional/circuits/bell_state.qasm 2

# Show the probability distribution and ASCII diagram
./qubby tests/functional/circuits/bell_state.qasm 2 --print --ascii
```

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

# Unit tests only
make unit_tests && ./tests/unit_tests

# Functional tests only
make functional_tests

# Generate HTML coverage report
make coverage
# Report available at coverage/html/index.html
```

## Benchmarking

Measures wall-clock time for `applyGate` (OpenMP) vs `applyGateSYCL` (SYCL) across 20–25 qubit counts.

```bash
# OpenMP benchmark
make benchmark

# SYCL benchmark (requires icpx)
make benchmark_sycl
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
├── src/                       # Implementation files
│   ├── main.cpp               # CLI entry point
│   ├── QuantumState.cpp       # State vector + gate logic
│   ├── QuantumState.sycl.cpp  # SYCL GPU kernel
│   ├── QuantumCircuit.cpp     # OpenQASM parser and circuit runner
│   ├── DeutschJozsa.cpp
│   ├── Grover.cpp
│   ├── QuantumFourierTransform.cpp
│   ├── Shor.cpp
│   ├── Simon.cpp
│   └── VariationalQuantumEigensolver.cpp
├── include/                   # Header files
├── tests/
│   ├── unit/                  # Criterion unit tests
│   └── functional/            # Shell-based functional tests + .qasm circuits
├── benchmarks/                # Benchmark sources and binaries
└── Makefile
```

## License

MIT
