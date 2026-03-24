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
#include <regex>

QuantumCircuit::QuantumCircuit(int nbQubits) {
    state.initialize(nbQubits);
}

/**
 * @details
 * Thin wrapper that delegates to QuantumState::measure().
 * The underlying state is collapsed after this call.
 */
int QuantumCircuit::measure() {
    return state.measure();
}

/**
 * @brief Extracts the integer qubit index from a token of the form `q[N]`.
 *
 * @param token String token to parse (e.g. `"q[2]"`).
 * @return The integer index N.
 * @throws std::runtime_error if the token does not contain matching `[` and `]`.
 */
int parseQubitIndexPos(std::string token) {
    size_t openPos = token.find('[');
    size_t closePos = token.find(']');
    if (openPos == std::string::npos || closePos == std::string::npos)
        throw std::runtime_error("Argument missing in q[...]");

    std::string content = token.substr(openPos+1, closePos-openPos-1);
    return std::stoi(content);
}

/**
 * @brief Extracts the rotation angle from a parametric gate token of the form `rx(θ)`.
 *
 * Supports three symbolic forms in addition to raw doubles:
 *   - `pi` → M_PI
 *   - `pi/N` → M_PI / N
 *   - `pi*N` → M_PI * N
 *
 * @param token Gate token including the angle (e.g. `"rx(pi/2)"`, `"rz(1.5708)"`).
 * @return The angle θ as a double (radians).
 * @throws std::runtime_error if the token does not contain matching `(` and `)`.
 */
double parseQubitIndexPar(std::string token) {
    size_t openPos = token.find('(');
    size_t closePos = token.find(')');
    if (openPos == std::string::npos || closePos == std::string::npos)
        throw std::runtime_error("Argument missing in r(...)");

    std::string content = token.substr(openPos+1, closePos-openPos-1);
    std::regex pattern_div(R"(pi\s*/\s*([0-9]+))");
    std::regex pattern_mult(R"(pi\s*\*\s*([0-9]+))");
    std::smatch matches;

    if (content == "pi") return M_PI;
    else if (std::regex_search(content, matches, pattern_div)) return M_PI / std::stoi(matches[1].str());
    else if (std::regex_search(content, matches, pattern_mult)) return M_PI * std::stoi(matches[1].str());
    return std::stod(content);
}

/**
 * @details
 * The parser reads the file line by line. For each line:
 * 1. Block-comment tracking: if `in_block_comment` is true, scan for `*\/` and skip
 *    everything up to (and including) it; otherwise scan for `/*` and strip the rest.
 * 2. Extract the first token as the gate name. Skip the line if it starts with `//`.
 * 3. For single-qubit gates (`h`, `x`, `z`): read one `q[N]` token, parse the index,
 *    and dispatch to the corresponding lambda in the `simple_gates` map.
 * 4. For parametric gates (`rx`, `ry`, `rz`): the angle θ is embedded in the gate
 *    token (`rx(pi/2)`), parsed by parseQubitIndexPar(), then dispatched via `multi_gates`.
 * 5. For `cx`: read two `q[N]` tokens (control, then target) and call cnotGate().
 * 6. For `swap`: read two `q[N]` tokens and call swapGate().
 * 7. For `ccx`: read three `q[N]` tokens (c0, c1, target) and call toffoliGate().
 * 8. For `qreg`: read one `q[N]` token and call initialize(N) to resize the state.
 * 9. Any other gate name throws std::runtime_error.
 */
void QuantumCircuit::load(std::string filename) {
    std::ifstream file(filename);

    if (!file.is_open())
        throw std::runtime_error("Error opening the file " + filename);

    std::string line;
    std::unordered_map<std::string, std::function<void(int)>> simple_gates = {
        {"h", [this](int q) { state.hGate(q); }},
        {"x", [this](int q) { state.xGate(q); }},
        {"z", [this](int q) { state.zGate(q); }},
    };
    std::unordered_map<std::string, std::function<void(double, int)>> multi_gates = {
        {"rx", [this](double t, int q) { state.rxGate(t, q); }},
        {"ry", [this](double t, int q) { state.ryGate(t, q); }},
        {"rz", [this](double t, int q) { state.rzGate(t, q); }},
    };

    bool in_block_comment = false;

    while (std::getline(file, line)) {
        if (in_block_comment) {
            size_t end_comment_pos = line.find("*/");
            if (end_comment_pos != std::string::npos) {
                line = line.substr(end_comment_pos + 2);
                in_block_comment = false;
            } else {
                continue;
            }
        }

        size_t start_comment_pos = line.find("/*");
        if (start_comment_pos != std::string::npos) {
            size_t end_comment_pos = line.find("*/", start_comment_pos + 2);
            if (end_comment_pos != std::string::npos) {
                line = line.substr(0, start_comment_pos) + line.substr(end_comment_pos + 2);
            } else {
                line = line.substr(0, start_comment_pos);
                in_block_comment = true;
            }
        }

        std::istringstream iss(line);
        std::string command;

        if (iss >> command && !command.starts_with("//")) {
            std::string token;
            iss >> token;
            int qubitIndex = parseQubitIndexPos(token);

            auto it_simple = simple_gates.find(command);
            std::string base = command.substr(0, command.find('('));
            auto it_multi = multi_gates.find(base);
            if (it_simple != simple_gates.end()) {
                it_simple->second(qubitIndex);
            } else if (it_multi != multi_gates.end()) {
                double theta = parseQubitIndexPar(command);
                it_multi->second(theta, qubitIndex);
            } else if (command == "cx") {
                std::string token2;
                iss >> token2;
                int qubitIndex2 = parseQubitIndexPos(token2);
                state.cnotGate(qubitIndex, qubitIndex2);
            } else if (command == "qreg") {
                int n = parseQubitIndexPos(token);
                state.initialize(n);
            } else {
                throw std::runtime_error("Gate " + command + " not found");
            }
        }
    }
}
