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
#include <algorithm>
#include <iostream>

QuantumCircuit::QuantumCircuit(int nbQubits) {
    state.initialize(nbQubits);
    this->numQubits = nbQubits;
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
                circuitLog.push_back({command, {qubitIndex}});
            } else if (it_multi != multi_gates.end()) {
                double theta = parseQubitIndexPar(command);
                it_multi->second(theta, qubitIndex);
                circuitLog.push_back({base, {qubitIndex}});
            } else if (command == "cx") {
                std::string token2;
                iss >> token2;
                int qubitIndex2 = parseQubitIndexPos(token2);
                state.cnotGate(qubitIndex, qubitIndex2);
                circuitLog.push_back({"cx", {qubitIndex, qubitIndex2}});
            } else if (command == "swap") {
                std::string token2;
                iss >> token2;
                int qubitIndex2 = parseQubitIndexPos(token2);
                state.swapGate(qubitIndex, qubitIndex2);
                circuitLog.push_back({"swap", {qubitIndex, qubitIndex2}});
            } else if (command == "ccx") {
                std::string token2, token3;
                iss >> token2 >> token3;
                int qubitIndex2 = parseQubitIndexPos(token2);
                int qubitIndex3 = parseQubitIndexPos(token3);
                state.toffoliGate(qubitIndex, qubitIndex2, qubitIndex3);
                circuitLog.push_back({"ccx", {qubitIndex, qubitIndex2, qubitIndex3}});
            } else if (command == "qreg") {
                int n = parseQubitIndexPos(token);
                state.initialize(n);
                numQubits = n;
                circuitLog.clear();
            } else {
                throw std::runtime_error("Gate " + command + " not found");
            }
        }
    }
}

/**
 * @details
 * Layout algorithm:
 * 1. Each qubit wire is a string initialised to `"q[i]: ---"`.
 * 2. For each gate in circuitLog:
 *    a. Pad all wires to the same column with `DASH` characters.
 *    b. Compute the column width as the maximum display width of any symbol
 *       that appears on a wire for this gate (minimum 1).
 *    c. For each wire: if it carries a symbol (gate operand or vertical connector),
 *       append `-<sym>-`; otherwise append `---` (pass-through dashes).
 * 3. After all gates, append three trailing dashes to every wire and print.
 *
 * The helper `vw()` counts Unicode code-points (not bytes) so that multi-byte
 * UTF-8 symbols (│, ●, ×) are treated as width 1 for alignment.
 * `innerSym()` maps (gate, qubit) → the symbol to draw, returning `""` for
 * uninvolved wires and `PIPE` for wires strictly between multi-qubit operands.
 */
void QuantumCircuit::draw() {
    static const std::string DASH  = "-";
    static const std::string PIPE  = "│";
    static const std::string DOT   = "●";
    static const std::string CROSS = "×";

    auto vw = [](const std::string& s) -> int {
        int n = 0;
        for (size_t i = 0; i < s.size(); ) {
            unsigned char c = (unsigned char)s[i];
            if (c < 0x80) i += 1;
            else if (c < 0xE0) i += 2;
            else if (c < 0xF0) i += 3;
            else i += 4;
            ++n;
        }
        return n;
    };

    auto rep = [](const std::string& ch, int n) {
        std::string r;
        for (int i = 0; i < n; ++i)
            r += ch;
        return r;
    };

    auto innerSym = [&](const GateRecord& g, int q) -> std::string {
        bool inGate = std::find(g.qubits.begin(), g.qubits.end(), q) != g.qubits.end();
        int  minQ = *std::min_element(g.qubits.begin(), g.qubits.end());
        int  maxQ = *std::max_element(g.qubits.begin(), g.qubits.end());
        bool between = !inGate && q > minQ && q < maxQ;

        if (!inGate && !between) return "";
        if (between) return PIPE;

        if (g.name == "h") return "H";
        if (g.name == "x") return "X";
        if (g.name == "z") return "Z";
        if (g.name == "rx") return "Rx";
        if (g.name == "ry") return "Ry";
        if (g.name == "rz") return "Rz";
        if (g.name == "cx") return (q == g.qubits[0]) ? DOT : "X";
        if (g.name == "swap") return CROSS;
        if (g.name == "ccx") return (q == g.qubits[2]) ? "X" : DOT;
        return "?";
    };

    std::vector<std::string> lines(numQubits);
    std::vector<int> vcols(numQubits, 0);
    for (int i = 0; i < numQubits; i++) {
        lines[i] = "q[" + std::to_string(i) + "]: " + rep(DASH, 3);
        vcols[i] = 3;
    }

    for (const auto& gate : circuitLog) {
        int maxCol = *std::max_element(vcols.begin(), vcols.end());
        for (int q = 0; q < numQubits; q++)
            while (vcols[q] < maxCol) {
                lines[q] += DASH;
                vcols[q]++;
            }

        int innerW = 1;
        for (int q = 0; q < numQubits; q++) {
            std::string sym = innerSym(gate, q);
            if (!sym.empty()) innerW = std::max(innerW, vw(sym));
        }

        for (int q = 0; q < numQubits; q++) {
            std::string sym = innerSym(gate, q);
            if (sym.empty()) {
                lines[q] += rep(DASH, innerW + 2);
            } else {
                int pad = innerW - vw(sym);
                lines[q] += DASH + sym + rep(DASH, pad + 1);
            }
            vcols[q] += innerW + 2;
        }
    }

    for (int q = 0; q < numQubits; q++)
        std::cout << lines[q] + rep(DASH, 3) << "\n";
}


void QuantumCircuit::printState() {
    state.printState();
}
