/*
 * Project: qubby
 * File name: main.cpp
 * Author: Cannelle Gourdet - lankley
 * File description: CLI entry point - loads and runs an OpenQASM circuit file,
 *                   then prints the measurement result.
 */

#include "QuantumCircuit.hpp"
#include "BlochSphere.hpp"
#include "REPL.hpp"
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <unordered_map>

int main(int argc, char *argv[]) {
    std::vector<std::string> positional;
    std::unordered_set<std::string> flags;
    std::unordered_map<std::string, std::string> kvflags;

    for (int i = 1; i < argc; i++) {
        std::string arg(argv[i]);
        if (arg.rfind("-", 0) == 0) {
            size_t eq = arg.find('=');
            if (eq != std::string::npos)
                kvflags[arg.substr(0, eq)] = arg.substr(eq + 1);
            else
                flags.insert(arg);
        } else
            positional.push_back(arg);
    }

    if (flags.count("-h") || flags.count("--help")) {
        std::cout << "Usage: " << argv[0] << " <circuit.qasm> <num_qubits> [OPTIONS]\n"
                  << "       " << argv[0] << " --repl [num_qubits]\n\n"
                  << "Arguments:\n"
                  << "  <circuit.qasm>           Path to an OpenQASM 2.0 circuit file\n"
                  << "  <num_qubits>             Number of qubits in the circuit\n\n"
                  << "Options:\n"
                  << "  --repl                   Start interactive OpenQASM REPL\n"
                  << "  --noise                  Enable realistic noise simulation (density matrix, ≤12 qubits)\n"
                  << "  --noise-rate=<p>         Per-gate depolarizing error rate in [0,1] (default: 0.01)\n"
                  << "  --ascii                  Draw an ASCII representation of the circuit\n"
                  << "  --print                  Print the quantum state before measurement\n"
                  << "  --bloch=<qubit>          Print Bloch sphere (ASCII) for the given qubit\n"
                  << "  --bloch-svg=<qubit>      Save Bloch sphere SVG (requires --bloch-file)\n"
                  << "  --bloch-file=<path>      Output file for Bloch sphere SVG\n"
                  << "  --export-qiskit=<file>   Export circuit as a Qiskit Python script\n"
                  << "  --export-cirq=<file>     Export circuit as a Cirq Python script\n"
                  << "  -h, --help               Show this help message and exit\n\n"
                  << "Noise channels (applied per qubit after each gate when --noise is active):\n"
                  << "  Depolarizing: ρ → (1-p)ρ + (p/3)(XρX + YρY + ZρZ)\n";
        return 0;
    }

    if (flags.count("--repl")) {
        int numQubits = 1;
        if (!positional.empty()) {
            try { numQubits = std::stoi(positional[0]); }
            catch (...) { numQubits = 1; }
        }
        try {
            REPL repl(numQubits);
            repl.run();
        } catch (const std::exception& e) {
            std::cerr << "Error: " << e.what() << "\n";
            return 1;
        }
        return 0;
    }

    if (positional.size() != 2) {
        std::cerr << "Usage: " << argv[0]
                  << " <circuit.qasm> <num_qubits> [OPTIONS]\n"
                  << "       " << argv[0] << " --repl [num_qubits]\n";
        return 1;
    }

    try {
        int numQubits = std::stoi(positional[1]);
        QuantumCircuit circuit(numQubits);

        if (flags.count("--noise")) {
            double rate = 0.01;
            if (kvflags.count("--noise-rate"))
                rate = std::stod(kvflags.at("--noise-rate"));
            circuit.validNoise(rate);
        }

        circuit.load(positional[0]);

        if (flags.count("--ascii"))
            circuit.draw();

        if (flags.count("--print"))
            circuit.printState();

        if (kvflags.count("--bloch") || kvflags.count("--bloch-svg")) {
            if (circuit.isNoisy()) {
                std::cerr << "Warning: Bloch sphere requires pure-state backend "
                             "(skip --noise to use --bloch).\n";
            } else {
                int blochQubit = 0;
                if (kvflags.count("--bloch"))
                    blochQubit = std::stoi(kvflags.at("--bloch"));
                else if (kvflags.count("--bloch-svg"))
                    blochQubit = std::stoi(kvflags.at("--bloch-svg"));

                auto bv = circuit.blochVector(blochQubit);
                BlochSphere::printASCII(bv, blochQubit);

                if (kvflags.count("--bloch-file"))
                    BlochSphere::saveSVG(bv, blochQubit, kvflags.at("--bloch-file"));
            }
        }

        if (kvflags.count("--export-qiskit"))
            circuit.exportQiskit(kvflags.at("--export-qiskit"));

        if (kvflags.count("--export-cirq"))
            circuit.exportCirq(kvflags.at("--export-cirq"));

        int result = circuit.measure();
        std::cout << "Measurement: " << result << "\n";

    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
