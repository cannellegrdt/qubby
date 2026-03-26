# qubby - blueprint

## 1. Technical architecture
* **Scale:** 20 to 25 Qubits (Memory: 16MB to 512MB RAM).
* **State representation:** A single contiguous `std::vector<std::complex<double>>` of size $2^n$.
* **Gate logic:** "In-place" transformations using **bitwise masking** (avoids massive matrix-matrix multiplications).
* **Acceleration:**
    * **CPU:** OpenMP (multithreading).
    * **GPU:** SYCL (via Intel oneAPI) to utilize my **Iris Xe** execution units.

---

## 2. The core "bitwise trick"
To apply a gate to qubit $k$ without a $2^n \times 2^n$ matrix:
1.  Iterate through all indices $i$ from $0$ to $2^{n-1}$.
2.  For each $i$, find the pair of states where the $k$-th bit is `0` (index $i_0$) and `1` (index $i_1$).
3.  Apply the $2 \times 2$ gate matrix to this pair:
    $$\begin{pmatrix} \psi[i_0] \\ \psi[i_1] \end{pmatrix} = \begin{pmatrix} U_{00} & U_{01} \\ U_{10} & U_{11} \end{pmatrix} \begin{pmatrix} \psi[i_0] \\ \psi[i_1] \end{pmatrix}$$

---

## 3. Detailed TODO list

### Phase 1: The foundation
* [x] **Project setup:** Configure Makefile with C++20 support.
* [x] **State vector class:** Create a `QuantumState` class that wraps `std::vector<std::complex<double>>`.
* [x] **Initialization:** Implement the `initialize(int num_qubits)` method (sets state to $|00...0\rangle$, where index 0 is $1.0$).
* [x] **Memory check:** Add a safety check to prevent allocating more RAM than your system has.

### Phase 2: Single-qubit gates
* [x] **Gate interface:** Create a generic `applyGate(int targetQubit, Matrix2x2 matrix)`.
* [x] **Standard gates:** Implement predefined functions for:
    * **X (NOT):** Flips the probability amplitudes.
    * **H (Hadamard):** Creates superposition.
    * **Z (Phase Flip):** Changes the sign of the $|1\rangle$ state.
* [x] **Validation:** Verify that after two H-gates, the state returns to $|0\rangle$.

### Phase 3: Entanglement & multi-qubit
* [x] **CNOT Implementation:** Implement the Controlled-NOT gate. This is the "soul" of quantum computing.
    * *Logic:* Only swap/transform indices where the *control* bit is `1`.
* [x] **Entanglement test:** Create a **bell state** $(\frac{|00\rangle + |11\rangle}{\sqrt{2}})$.
* [x] **Measurement:** Implement a `measure()` function using `std::discrete_distribution` to collapse the state based on $| \alpha |^2$.

### Phase 4: OpenQASM parser
* [x] **Lexer:** Write a simple string parser to read `.qasm` files line by line.
* [x] **Command mapping:** Map strings like `"h q[0]"` to your C++ function `state.hGate(0)`.
* [x] **Circuit runner:** Create a `QuantumCircuit` class that loads a file and executes the gates in sequence.

### Phase 5: Hardware acceleration
* [x] **OpenMP (CPU):** Add `#pragma omp parallel for` to your gate application loops. This will immediately use all your CPU threads.
* [x] **SYCL (GPU):** Install **Intel oneAPI Base Toolkit**.
* [x] **Compute kernel:** Write the gate application loop as a SYCL kernel to offload the work to the Iris Xe's 96 Execution Units.

### Phase 6: Parametric gates & extended gate set
* [x] **Rotation gates:** Implement `Rx(θ)`, `Ry(θ)`, `Rz(θ)` — arbitrary single-qubit rotations around the Bloch sphere axes.
* [x] **Multi-qubit gates:** Implement `swapGate(int q0, int q1)` and `toffoliGate(int c0, int c1, int target)` (Controlled-Controlled-NOT).
* [x] **OpenQASM 2.0 support:** Extend the parser to handle register declarations (`qreg q[n];`), parametric gates (`rx(pi/2) q[0]`), and multi-line comments.

### Phase 7: Quantum algorithms
* [x] **Deutsch-Jozsa:** Determine whether a function is constant or balanced with a single query — the first demonstration of quantum advantage.
* [x] **Grover's algorithm:** Implement amplitude amplification to search an unstructured database in $O(\sqrt{N})$ instead of $O(N)$.
* [x] **Quantum Fourier Transform (QFT):** Implement QFT as a building block for more advanced algorithms (e.g. Shor's).

### Phase 8: Observability & benchmarking
* [x] **Probability distribution display:** Print the amplitude distribution of all basis states after a circuit run.
* [x] **ASCII circuit drawer:** Render the gate sequence as a readable ASCII diagram in the terminal.
* [x] **OpenMP vs SYCL benchmark:** Measure and compare wall-clock time for `applyGate` vs `applyGateSYCL` on large qubit counts (20–25 qubits).

### Phase 9: Advanced quantum algorithms
* [x] **Shor's algorithm:** Implement integer factorisation in polynomial time — uses QFT as its core subroutine.
* [x] **Simon's algorithm:** Find the period of a function with a single oracle call — exponentially faster than any classical algorithm.
* [x] **Variational Quantum Eigensolver (VQE):** Hybrid classical/quantum algorithm to find ground-state energies using parametric circuits.

### Phase 10: Realistic noise simulation
* [ ] **Noise channels:** Model quantum decoherence with bit-flip, phase-flip, and depolarizing channels.
* [ ] **Density matrix:** Replace the state vector with a density matrix to simulate mixed states.
* [ ] **Imperfect gates:** Add an error parameter to gates to simulate real hardware noise.

### Phase 11: Visualisation & interoperability
* [ ] **Bloch sphere:** Render a single-qubit state on the Bloch sphere (SVG or terminal graphics).
* [ ] **Export to Qiskit/Cirq:** Generate equivalent Python code from a loaded circuit for cross-platform validation.
* [ ] **Interactive REPL:** A shell that accepts OpenQASM instructions in real time and updates the state live.
