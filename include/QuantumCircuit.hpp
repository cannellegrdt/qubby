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

        /** @brief Number of qubits in the circuit (updated by `qreg`). */
        int numQubits;

        /**
         * @struct GateRecord
         * @brief Stores one gate application for later rendering by draw().
         *
         * @var name   Gate name as it appears in the OpenQASM source (e.g. `"h"`, `"cx"`).
         * @var qubits Ordered list of qubit indices the gate acts on.
         *             Single-qubit gates have one entry; cx/swap have two; ccx has three.
         */
        struct GateRecord {
            std::string name;
            std::vector<int> qubits;
        };

        /** @brief Sequential log of every gate applied since the last `qreg` (or construction). */
        std::vector<GateRecord> circuitLog;

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

        /**
         * @brief Renders the circuit as an ASCII diagram to stdout.
         *
         * Replays the gate log (`circuitLog`) and produces one horizontal wire per qubit.
         * Gate symbols are aligned in columns; multi-qubit gates draw a vertical connector
         * (`│`) through any qubits between the operands.
         *
         * Symbol legend:
         * @code
         *   H  Ry  Rz  Rx      single-qubit gates
         *   ●                  control qubit (cx, ccx)
         *   X                  target qubit (cx, ccx)
         *   ×                  SWAP operand
         *   │                  vertical wire between multi-qubit gate operands
         * @endcode
         *
         * Unicode characters (│, ●, ×) are counted by display width, not byte length,
         * so columns stay aligned even with multi-byte UTF-8 symbols.
         */
        void draw();

        /**
         * @brief Prints the non-negligible amplitudes of the current quantum state to stdout.
         *
         * Delegates to QuantumState::printState(). Only basis states with
         * probability > 1e-6 are shown.
         */
        void printState();
};

#endif /* QUANTUMCIRCUIT_HPP_ */
