/*
 * Project: qubby
 * File name: QuantumCircuit.cpp
 * Author: Cannelle Gourdet - lankley
 * File description: Implements the QuantumCircuit methods: constructor and OpenQASM file parser.
 */

#include "QuantumCircuit.hpp"
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <unordered_map>
#include <functional>

QuantumCircuit::QuantumCircuit(int nbQubits) {
    state.initialize(nbQubits);
}

/**
 * @brief Extracts the integer qubit index from a token of the form `q[N]`.
 *
 * @param token String token to parse (e.g. `"q[2]"`).
 * @return The integer index N.
 * @throws std::runtime_error if the token does not contain matching `[` and `]`.
 */
int parseQubitIndex(std::string token) {
    size_t openPos = token.find('[');
    size_t closePos = token.find(']');
    if (openPos == std::string::npos || closePos == std::string::npos)
        throw std::runtime_error("Argument missing in q[...]");

    std::string content = token.substr(openPos+1, closePos-openPos-1);
    return std::stoi(content);
}

/**
 * @details
 * The parser reads the file line by line. For each line:
 * 1. Extract the first token as the gate name. Skip the line if it starts with `//`.
 * 2. For single-qubit gates (`h`, `x`, `z`): read one `q[N]` token, parse the index,
 *    and dispatch to the corresponding lambda in the `gates` map.
 * 3. For `cx`: read two `q[N]` tokens (control, then target) and call cnotGate().
 * 4. Any other gate name throws std::runtime_error.
 *
 * The gate dispatch table is built as an unordered_map of lambdas so that adding
 * new single-qubit gates only requires inserting one entry.
 */
void QuantumCircuit::load(std::string filename) {
    std::ifstream file(filename);

    if (!file.is_open())
        throw std::runtime_error("Error opening the file " + filename);

    std::string line;
    std::unordered_map<std::string, std::function<void(int)>> gates = {
                {"h", [this](int q) { state.hGate(q); }},
                {"x", [this](int q) { state.xGate(q); }},
                {"z", [this](int q) { state.zGate(q); }},
            };

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string command;

        if (iss >> command && !command.starts_with("//")) {
            std::string token;
            iss >> token;
            int qubitIndex = parseQubitIndex(token);

            
            auto it = gates.find(command);
            if (it != gates.end()) {
                it->second(qubitIndex);
            } else if (command == "cx") {
                std::string token2;
                iss >> token2;
                int qubitIndex2 = parseQubitIndex(token2);
                state.cnotGate(qubitIndex, qubitIndex2);
            } else {
                throw std::runtime_error("Gate " + command + " not found");
            }
        }
    }
}
