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

int main(int argc, char *argv[]) {
    std::vector<std::string> positional;
    std::unordered_set<std::string> flags;

    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        if (arg.rfind("-", 0) == 0)
            flags.insert(arg);
        else
            positional.push_back(arg);
    }

    if (flags.count("-h") || flags.count("--help")) {
        std::cout << "Usage: " << argv[0] << " <circuit.qasm> <num_qubits> [OPTIONS]\n\n"
                  << "Arguments:\n"
                  << "  <circuit.qasm>   Path to an OpenQASM 2.0 circuit file\n"
                  << "  <num_qubits>     Number of qubits in the circuit\n\n"
                  << "Options:\n"
                  << "  --ascii          Draw an ASCII representation of the circuit\n"
                  << "  --print          Print the quantum state before measurement\n"
                  << "  -h, --help       Show this help message and exit\n";
        return 0;
    }

    if (positional.size() != 2) {
        std::cerr << "Usage: " << argv[0] << " <circuit.qasm> <num_qubits> [--ascii] [--print]\n";
        return 1;
    }

    try {
        int numQubits = std::stoi(positional[1]);
        QuantumCircuit circuit(numQubits);
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
