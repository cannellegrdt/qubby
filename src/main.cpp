/*
 * Project: qubby
 * File name: main.cpp
 * Author: Cannelle Gourdet - lankley
 * File description: CLI entry point - loads and runs an OpenQASM circuit file,
 *                   then prints the measurement result.
 */

#include "QuantumCircuit.hpp"
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
        std::cout << "Usage: " << argv[0] << " <circuit.qasm> <num_qubits> [OPTIONS]\n\n"
                  << "Arguments:\n"
                  << "  <circuit.qasm>        Path to an OpenQASM 2.0 circuit file\n"
                  << "  <num_qubits>          Number of qubits in the circuit\n\n"
                  << "Options:\n"
                  << "  --noise               Enable realistic noise simulation (density matrix, ≤12 qubits)\n"
                  << "  --noise-rate=<p>      Per-gate depolarizing error rate in [0,1] (default: 0.01)\n"
                  << "  --ascii               Draw an ASCII representation of the circuit\n"
                  << "  --print               Print the quantum state before measurement\n"
                  << "  -h, --help            Show this help message and exit\n\n"
                  << "Noise channels (applied per qubit after each gate when --noise is active):\n"
                  << "  Depolarizing: ρ → (1-p)ρ + (p/3)(XρX + YρY + ZρZ)\n";
        return 0;
    }

    if (positional.size() != 2) {
        std::cerr << "Usage: " << argv[0] << " <circuit.qasm> <num_qubits> [--noise] [--noise-rate=<p>] [--ascii] [--print]\n";
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

        int result = circuit.measure();
        std::cout << "Measurement: " << result << "\n";
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
