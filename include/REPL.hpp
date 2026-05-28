/*
 * Project: qubby
 * File name: REPL.hpp
 * Author: Cannelle Gourdet - lankley
 * File description: Declares the REPL class — an interactive shell that accepts OpenQASM
 *                   instructions in real time and updates the quantum state live.
 */

#ifndef REPL_HPP_
    #define REPL_HPP_

    #include "QuantumCircuit.hpp"
    #include "BlochSphere.hpp"
    #include <string>

/**
 * @class REPL
 * @brief Interactive OpenQASM Read-Eval-Print Loop.
 *
 * Accepts one OpenQASM instruction per line. After each instruction the quantum
 * state is printed automatically. Special dot-commands control the session:
 *
 *   .help                    — show this list
 *   .quit  /  .exit          — exit the REPL
 *   .state                   — print the current probability distribution
 *   .measure                 — collapse and print measurement result
 *   .reset                   — reinitialise to |0...0⟩ (same qubit count)
 *   .ascii                   — draw the ASCII circuit diagram
 *   .bloch <qubit>           — render Bloch sphere in terminal (ASCII)
 *   .bloch-svg <qubit> <f>   — save Bloch sphere SVG to file f
 *   .export-qiskit <file>    — export circuit to Qiskit Python script
 *   .export-cirq <file>      — export circuit to Cirq Python script
 *   .noise <rate>            — enable depolarizing noise (rate in [0,1])
 *
 * Regular OpenQASM lines (e.g. `h q[0]`, `cx q[0] q[1]`, `qreg q[3]`) are
 * dispatched directly to QuantumCircuit::executeInstruction().
 */
class REPL {
public:
    /**
     * @brief Constructs the REPL, initialised with @p numQubits qubits.
     * @param numQubits Starting qubit count (default 1). The user can change
     *                  this at any time via `qreg q[N]`.
     */
    explicit REPL(int numQubits = 1);

    /**
     * @brief Starts the interactive loop; returns when the user types `.quit`.
     *
     * Reads lines from stdin. EOF (Ctrl-D) is treated as `.quit`.
     */
    void run();

private:
    QuantumCircuit circuit;

    /**
     * @brief Processes one input line (QASM instruction or dot-command).
     * @return false if the session should end (`.quit` / `.exit`).
     */
    bool processLine(const std::string& line);

    void printHelp() const;
    void printPrompt() const;
};

#endif /* REPL_HPP_ */
