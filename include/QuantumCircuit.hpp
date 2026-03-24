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
 * Lines starting with `//` are treated as comments and ignored.
 *
 * Supported gates: `h`, `x`, `z`, `cx`.
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
         * @brief Parses and executes a circuit file line by line.
         *
         * Each non-comment line must follow the format:
         * @code
         *   <gate> q[<index>]              // single-qubit gates: h, x, z
         *   cx q[<control>] q[<target>]    // two-qubit gate
         * @endcode
         *
         * @param filename Path to the circuit file.
         * @throws std::runtime_error if the file cannot be opened, a qubit index
         *         is malformed, or an unknown gate name is encountered.
         */
        void load(std::string filename);
};

#endif /* QUANTUMCIRCUIT_HPP_ */
