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

int main(int argc, char *argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <circuit.qasm> <num_qubits>\n";
        return 1;
    }

    try {
        int numQubits = std::stoi(argv[2]);
        QuantumCircuit circuit(numQubits);
        circuit.load(argv[1]);
        int result = circuit.measure();
        std::cout << "Measurement: " << result << "\n";
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
