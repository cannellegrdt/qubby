/*
 * Project: qubby
 * File name: QuantumCircuit.hpp
 * Author: Cannelle Gourdet - lankley
 * File description: Declares the QuantumCircuit class, which loads and executes a quantum circuit
 *                   from an OpenQASM file, dispatching each instruction to a QuantumState.
 */

#ifndef QUANTUMCIRCUIT_HPP_
    #define QUANTUMCIRCUIT_HPP_

    #include "QuantumState.hpp"
    #include <string>

/**
 * @class QuantumCircuit
 * @brief Loads and runs a quantum circuit described in an OpenQASM text file.
 *
 * Wraps a QuantumState and exposes a simple file-based interface. Each line of
 * the input file is parsed as a gate instruction (e.g. `h q[0]`, `cx q[0] q[1]`).
 *
 * Comment styles supported:
 *   - Line comments: `// ...`
 *   - Block comments: `/* ... *\/` (spanning multiple lines)
 *
 * Supported instructions:
 *   - `qreg q[N]`                 reinitialise state to N qubits
 *   - `h|x|z q[N]`                single-qubit gates
 *   - `cx q[control] q[target]`   CNOT
 *   - `swap q[q0] q[q1]`          SWAP
 *   - `ccx q[c0] q[c1] q[target]` Toffoli
 *   - `rx(θ)|ry(θ)|rz(θ) q[N]`    rotation gates (θ in radians, or symbolic: `pi`, `pi/N`, `pi*N`)
 */
class QuantumCircuit {
    private:
        /** @brief The underlying quantum state that gates are applied to. */
        QuantumState state;

    public:
        /**
         * @brief Constructs a QuantumCircuit and initialises the state to |00...0⟩.
         *
         * @param nbQubits Number of qubits in the circuit.
         * @throws std::runtime_error if the required memory exceeds MB_LIMIT.
         */
        QuantumCircuit(int nbQubits);

        /**
         * @brief Measures the underlying quantum state and collapses it.
         *
         * @return The measured basis-state index in [0, 2^n).
         */
        int measure();

        /**
         * @brief Parses and executes a circuit file line by line.
         *
         * Accepts a subset of OpenQASM 2.0 syntax:
         * @code
         *   qreg q[N];                         // (re)initialise to N qubits
         *   h|x|z q[N];                        // single-qubit gates
         *   cx q[control] q[target];           // CNOT
         *   swap q[q0] q[q1];                  // SWAP
         *   ccx q[c0] q[c1] q[target];         // Toffoli
         *   rx(pi/2)|ry(pi)|rz(1.5708) q[N];   // rotation gates
         *   // line comment
         *   /* block comment *\/
         * @endcode
         *
         * @param filename Path to the `.qasm` circuit file.
         * @throws std::runtime_error if the file cannot be opened, a qubit index
         *         is malformed, or an unknown gate name is encountered.
         */
        void load(std::string filename);
};

#endif /* QUANTUMCIRCUIT_HPP_ */
