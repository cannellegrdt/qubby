/*
 * Project: qubby
 * File name: REPL.cpp
 * Author: Cannelle Gourdet - lankley
 * File description: Implements the interactive OpenQASM REPL.
 */

#include "REPL.hpp"
#include <iostream>
#include <sstream>
#include <string>
#include <stdexcept>

REPL::REPL(int numQubits) : circuit(numQubits) {}

void REPL::printPrompt() const {
    std::cout << "qubby> " << std::flush;
}

void REPL::printHelp() const {
    std::cout <<
        "\n  OpenQASM instructions (same syntax as .qasm files):\n"
        "    qreg q[N]              — (re)initialise to N qubits\n"
        "    h|x|z q[N]             — single-qubit gate\n"
        "    rx(θ)|ry(θ)|rz(θ) q[N] — rotation gate (θ in radians or pi/N)\n"
        "    cx q[c] q[t]           — CNOT\n"
        "    swap q[a] q[b]         — SWAP\n"
        "    ccx q[c0] q[c1] q[t]   — Toffoli\n"
        "\n  Session commands:\n"
        "    .help                  — show this help\n"
        "    .quit  /  .exit        — exit the REPL\n"
        "    .state                 — print probability distribution\n"
        "    .measure               — collapse state, print result\n"
        "    .reset                 — reinitialise to |0...0>\n"
        "    .ascii                 — draw ASCII circuit diagram\n"
        "    .bloch <q>             — Bloch sphere terminal view for qubit q\n"
        "    .bloch-svg <q> <file>  — save Bloch sphere SVG for qubit q\n"
        "    .export-qiskit <file>  — export circuit as Qiskit Python script\n"
        "    .export-cirq <file>    — export circuit as Cirq Python script\n"
        "    .noise <rate>          — enable depolarizing noise (rate in [0,1])\n"
        "\n";
}

/**
 * @details
 * Dot-commands are dispatched by prefix matching on the first token.
 * Regular lines are forwarded to QuantumCircuit::executeInstruction().
 * After each successful gate instruction the state is printed automatically
 * (mirroring the live-update behaviour expected of a REPL).
 * Errors are caught and reported without terminating the session.
 */
bool REPL::processLine(const std::string& raw) {
    std::string line = raw;
    while (!line.empty() && (line.front() == ' ' || line.front() == '\t'))
        line.erase(line.begin());
    while (!line.empty() && (line.back() == ' ' || line.back() == '\t' || line.back() == '\r'))
        line.pop_back();

    if (line.empty()) return true;

    if (line[0] == '.') {
        std::istringstream iss(line);
        std::string cmd;
        iss >> cmd;

        if (cmd == ".quit" || cmd == ".exit") {
            std::cout << "Bye.\n";
            return false;
        }

        if (cmd == ".help") {
            printHelp();
            return true;
        }

        if (cmd == ".state") {
            circuit.printState();
            return true;
        }

        if (cmd == ".measure") {
            try {
                int result = circuit.measure();
                std::cout << "  Measured: " << result << " (|";
                int n = circuit.getNumQubits();
                for (int b = n - 1; b >= 0; b--)
                    std::cout << ((result >> b) & 1);
                std::cout << "\xe2\x9f\xa9" << ")\n";  // ⟩
            } catch (const std::exception& e) {
                std::cerr << "  Error: " << e.what() << "\n";
            }
            return true;
        }

        if (cmd == ".reset") {
            try {
                int n = circuit.getNumQubits();
                circuit.executeInstruction("qreg q[" + std::to_string(n) + "]");
                std::cout << "  Reset to |0...0\xe2\x9f\xa9 (" << n << " qubits)\n";
            } catch (const std::exception& e) {
                std::cerr << "  Error: " << e.what() << "\n";
            }
            return true;
        }

        if (cmd == ".ascii") {
            circuit.draw();
            return true;
        }

        if (cmd == ".bloch") {
            int q = 0;
            if (!(iss >> q)) {
                std::cerr << "  Usage: .bloch <qubit>\n";
                return true;
            }
            try {
                if (circuit.isNoisy()) {
                    std::cerr << "  Bloch sphere requires pure-state backend (no .noise)\n";
                    return true;
                }
                auto bv = circuit.blochVector(q);
                BlochSphere::printASCII(bv, q);
            } catch (const std::exception& e) {
                std::cerr << "  Error: " << e.what() << "\n";
            }
            return true;
        }

        if (cmd == ".bloch-svg") {
            int q = 0;
            std::string svgfile;
            if (!(iss >> q >> svgfile)) {
                std::cerr << "  Usage: .bloch-svg <qubit> <file.svg>\n";
                return true;
            }
            try {
                if (circuit.isNoisy()) {
                    std::cerr << "  Bloch sphere requires pure-state backend (no .noise)\n";
                    return true;
                }
                auto bv = circuit.blochVector(q);
                BlochSphere::printASCII(bv, q);
                BlochSphere::saveSVG(bv, q, svgfile);
            } catch (const std::exception& e) {
                std::cerr << "  Error: " << e.what() << "\n";
            }
            return true;
        }

        if (cmd == ".export-qiskit") {
            std::string fname;
            if (!(iss >> fname)) {
                std::cerr << "  Usage: .export-qiskit <file.py>\n";
                return true;
            }
            try {
                circuit.exportQiskit(fname);
            } catch (const std::exception& e) {
                std::cerr << "  Error: " << e.what() << "\n";
            }
            return true;
        }

        if (cmd == ".export-cirq") {
            std::string fname;
            if (!(iss >> fname)) {
                std::cerr << "  Usage: .export-cirq <file.py>\n";
                return true;
            }
            try {
                circuit.exportCirq(fname);
            } catch (const std::exception& e) {
                std::cerr << "  Error: " << e.what() << "\n";
            }
            return true;
        }

        if (cmd == ".noise") {
            double rate = 0.01;
            iss >> rate;
            try {
                circuit.validNoise(rate);
                std::cout << "  Noise enabled (depolarizing rate = " << rate << ")\n";
            } catch (const std::exception& e) {
                std::cerr << "  Error: " << e.what() << "\n";
            }
            return true;
        }

        std::cerr << "  Unknown command: " << cmd << "  (type .help)\n";
        return true;
    }

    try {
        circuit.executeInstruction(line);
        circuit.printState();
    } catch (const std::exception& e) {
        std::cerr << "  Error: " << e.what() << "\n";
    }
    return true;
}

/**
 * @details
 * Reads lines from stdin until EOF or a `.quit` / `.exit` command.
 * Prints a welcome banner and the initial state before the first prompt.
 */
void REPL::run() {
    std::cout << "\nqubby REPL  —  type OpenQASM instructions, .help for help, .quit to exit\n";
    std::cout << "Initial state: " << circuit.getNumQubits() << " qubit(s) in |0...0\xe2\x9f\xa9\n\n";

    std::string line;
    while (true) {
        printPrompt();
        if (!std::getline(std::cin, line)) {
            std::cout << "\nBye.\n";
            break;
        }
        if (!processLine(line))
            break;
    }
}
