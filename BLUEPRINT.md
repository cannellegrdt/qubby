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
