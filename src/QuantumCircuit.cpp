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
#include <cmath>

QuantumCircuit::QuantumCircuit(int nbQubits) {
    state.initialize(nbQubits);
    this->numQubits = nbQubits;
}

/**
 * @details
 * Thin wrapper that delegates to QuantumState::measure().
 * The underlying state is collapsed after this call.
 */
/**
 * @details
 * When noise is active, delegates to DensityMatrix::measure() which samples
 * from the diagonal of ρ (Born rule for mixed states) and collapses to |result⟩⟨result|.
 * When noise is inactive, delegates to QuantumState::measure() as before.
 */
int QuantumCircuit::measure() {
    if (noiseActivated)
        return dm.measure();
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
 *    everything up to (and including) it; otherwise scan for `/*` and strip the rest. // NOLINT
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

 * When noiseActivated is true, every gate is dispatched to the DensityMatrix backend
 * (dm) instead of QuantumState, and a depolarizing channel is applied to each affected
 * qubit immediately after the gate. This models imperfect hardware where each physical
 * gate introduces a small error with probability noiseModel.depolarizeRate.
 *
 * The circuitLog is always updated regardless of the noise mode, so draw() works in
 * both pure and noisy modes.
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
    std::unordered_map<std::string, std::function<void(int)>> dm_simple_gates = {
        {"h", [this](int q) { dm.hGate(q); }},
        {"x", [this](int q) { dm.xGate(q); }},
        {"z", [this](int q) { dm.zGate(q); }},
    };
    std::unordered_map<std::string, std::function<void(double, int)>> multi_gates = {
        {"rx", [this](double t, int q) { state.rxGate(t, q); }},
        {"ry", [this](double t, int q) { state.ryGate(t, q); }},
        {"rz", [this](double t, int q) { state.rzGate(t, q); }},
    };
    std::unordered_map<std::string, std::function<void(double, int)>> dm_multi_gates = {
        {"rx", [this](double t, int q) { dm.rxGate(t, q); }},
        {"ry", [this](double t, int q) { dm.ryGate(t, q); }},
        {"rz", [this](double t, int q) { dm.rzGate(t, q); }},
    };

    auto applyNoise = [this](std::vector<int> qubits) {
        if (!noiseActivated || noiseModel.depolarizeRate <= 0.0)
            return;
        for (int q : qubits)
            dm.applyDepolarizing(q, noiseModel.depolarizeRate);
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
                if (noiseActivated) {
                    dm_simple_gates[command](qubitIndex);
                    applyNoise({qubitIndex});
                } else
                    it_simple->second(qubitIndex);
                circuitLog.push_back({command, {qubitIndex}});

            } else if (it_multi != multi_gates.end()) {
                double theta = parseQubitIndexPar(command);
                if (noiseActivated) {
                    dm_multi_gates[base](theta, qubitIndex);
                    applyNoise({qubitIndex});
                } else
                    it_multi->second(theta, qubitIndex);
                circuitLog.push_back({base, {qubitIndex}, theta});

            } else if (command == "cx") {
                std::string token2;
                iss >> token2;
                int qubitIndex2 = parseQubitIndexPos(token2);
                if (noiseActivated) {
                    dm.cnotGate(qubitIndex, qubitIndex2);
                    applyNoise({qubitIndex, qubitIndex2});
                } else
                    state.cnotGate(qubitIndex, qubitIndex2);
                circuitLog.push_back({"cx", {qubitIndex, qubitIndex2}});

            } else if (command == "swap") {
                std::string token2;
                iss >> token2;
                int qubitIndex2 = parseQubitIndexPos(token2);
                if (noiseActivated) {
                    dm.swapGate(qubitIndex, qubitIndex2);
                    applyNoise({qubitIndex, qubitIndex2});
                } else
                    state.swapGate(qubitIndex, qubitIndex2);
                circuitLog.push_back({"swap", {qubitIndex, qubitIndex2}});

            } else if (command == "ccx") {
                std::string token2, token3;
                iss >> token2 >> token3;
                int qubitIndex2 = parseQubitIndexPos(token2);
                int qubitIndex3 = parseQubitIndexPos(token3);
                if (noiseActivated) {
                    dm.toffoliGate(qubitIndex, qubitIndex2, qubitIndex3);
                    applyNoise({qubitIndex, qubitIndex2, qubitIndex3});
                } else
                    state.toffoliGate(qubitIndex, qubitIndex2, qubitIndex3);
                circuitLog.push_back({"ccx", {qubitIndex, qubitIndex2, qubitIndex3}});

            } else if (command == "qreg") {
                int n = parseQubitIndexPos(token);
                state.initialize(n);
                numQubits = n;
                circuitLog.clear();
                if (noiseActivated) {
                    if (n > DM_MAX_QUBITS)
                        throw std::runtime_error(
                            "Noise simulation limited to " + std::to_string(DM_MAX_QUBITS) + " qubits");
                    dm.initialize(n);
                }
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
    static const std::string DASH = "-";
    static const std::string PIPE = "│";
    static const std::string DOT = "●";
    static const std::string CROSS = "×";

    auto vw = [](const std::string& s) -> int {
        int n = 0;
        for (size_t i = 0; i < s.size(); ) {
            unsigned char c = (unsigned char)s[i];
            if (c < 0x80) i += 1;
            else if (c < 0xE0) i += 2;
            else if (c < 0xF0) i += 3;
            else i += 4;
            n++;
        }
        return n;
    };

    auto rep = [](const std::string& ch, int n) {
        std::string r;
        for (int i = 0; i < n; i++)
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
    if (noiseActivated) {
        dm.printState();
        return;
    }
    state.printState();
}

/**
 * @details
 * Initialises the DensityMatrix backend to |0...0⟩⟨0...0| and activates the
 * noise flag so that load() routes all subsequent gate applications through dm
 * and appends a per-qubit depolarizing channel after each gate.
 *
 * Must be called before load() — gate applications that occurred before this call
 * will not be reflected in the density matrix.
 */
int QuantumCircuit::getNumQubits() const { return numQubits; }
bool QuantumCircuit::isNoisy() const { return noiseActivated; }

BlochSphere::BlochVector QuantumCircuit::blochVector(int qubit) const {
    if (noiseActivated)
        throw std::runtime_error(
            "Bloch sphere requires pure-state backend — do not use --noise");
    return BlochSphere::compute(state, qubit);
}

/**
 * @details
 * Mirrors the parsing logic of load() for a single line. Block-comment stripping
 * is intentionally omitted (the REPL processes one line at a time). Steps:
 *  1. Strip trailing `;`, line-comment `//...`, and surrounding whitespace.
 *  2. Extract the first token as the gate name.
 *  3. Dispatch to the same gate/DM backends as load(), and update circuitLog.
 */
void QuantumCircuit::executeInstruction(const std::string& line) {
    std::string stripped = line;

    auto cpos = stripped.find("//");
    if (cpos != std::string::npos)
        stripped = stripped.substr(0, cpos);

    while (!stripped.empty() &&
           (stripped.back() == ';' || stripped.back() == ' ' || stripped.back() == '\t'))
        stripped.pop_back();

    if (stripped.empty()) return;

    std::istringstream iss(stripped);
    std::string command;
    if (!(iss >> command) || command.starts_with("//")) return;

    std::string token;
    iss >> token;

    auto applyNoise = [this](std::vector<int> qubits) {
        if (!noiseActivated || noiseModel.depolarizeRate <= 0.0) return;
        for (int q : qubits)
            dm.applyDepolarizing(q, noiseModel.depolarizeRate);
    };

    std::string base = command.substr(0, command.find('('));

    if (command == "h" || command == "x" || command == "z") {
        int q = parseQubitIndexPos(token);
        if (noiseActivated) {
            if (command == "h") dm.hGate(q);
            else if (command == "x") dm.xGate(q);
            else dm.zGate(q);
            applyNoise({q});
        } else {
            if (command == "h") state.hGate(q);
            else if (command == "x") state.xGate(q);
            else state.zGate(q);
        }
        circuitLog.push_back({command, {q}});

    } else if (base == "rx" || base == "ry" || base == "rz") {
        int q = parseQubitIndexPos(token);
        double theta = parseQubitIndexPar(command);
        if (noiseActivated) {
            if (base == "rx") dm.rxGate(theta, q);
            else if (base == "ry") dm.ryGate(theta, q);
            else dm.rzGate(theta, q);
            applyNoise({q});
        } else {
            if (base == "rx") state.rxGate(theta, q);
            else if (base == "ry") state.ryGate(theta, q);
            else state.rzGate(theta, q);
        }
        circuitLog.push_back({base, {q}, theta});

    } else if (command == "cx") {
        std::string token2;
        iss >> token2;
        int q0 = parseQubitIndexPos(token);
        int q1 = parseQubitIndexPos(token2);
        if (noiseActivated) { dm.cnotGate(q0, q1); applyNoise({q0, q1}); }
        else state.cnotGate(q0, q1);
        circuitLog.push_back({"cx", {q0, q1}});

    } else if (command == "swap") {
        std::string token2;
        iss >> token2;
        int q0 = parseQubitIndexPos(token);
        int q1 = parseQubitIndexPos(token2);
        if (noiseActivated) { dm.swapGate(q0, q1); applyNoise({q0, q1}); }
        else state.swapGate(q0, q1);
        circuitLog.push_back({"swap", {q0, q1}});

    } else if (command == "ccx") {
        std::string token2, token3;
        iss >> token2 >> token3;
        int q0 = parseQubitIndexPos(token);
        int q1 = parseQubitIndexPos(token2);
        int q2 = parseQubitIndexPos(token3);
        if (noiseActivated) { dm.toffoliGate(q0, q1, q2); applyNoise({q0, q1, q2}); }
        else state.toffoliGate(q0, q1, q2);
        circuitLog.push_back({"ccx", {q0, q1, q2}});

    } else if (command == "qreg") {
        int n = parseQubitIndexPos(token);
        state.initialize(n);
        numQubits = n;
        circuitLog.clear();
        if (noiseActivated) {
            if (n > DM_MAX_QUBITS)
                throw std::runtime_error(
                    "Noise simulation limited to " + std::to_string(DM_MAX_QUBITS) + " qubits");
            dm.initialize(n);
        }
    } else {
        throw std::runtime_error("Gate " + command + " not found");
    }
}

void QuantumCircuit::exportQiskit(const std::string& filename) const {
    std::ofstream f(filename);
    if (!f.is_open())
        throw std::runtime_error("Cannot open file: " + filename);

    f << "from qiskit import QuantumCircuit\n\n";
    f << "qc = QuantumCircuit(" << numQubits << ")\n";

    for (const auto& g : circuitLog) {
        if (g.name == "h")
            f << "qc.h("    << g.qubits[0] << ")\n";
        else if (g.name == "x")
            f << "qc.x("    << g.qubits[0] << ")\n";
        else if (g.name == "z")
            f << "qc.z("    << g.qubits[0] << ")\n";
        else if (g.name == "rx")
            f << "qc.rx("   << g.param << ", " << g.qubits[0] << ")\n";
        else if (g.name == "ry")
            f << "qc.ry("   << g.param << ", " << g.qubits[0] << ")\n";
        else if (g.name == "rz")
            f << "qc.rz("   << g.param << ", " << g.qubits[0] << ")\n";
        else if (g.name == "cx")
            f << "qc.cx("   << g.qubits[0] << ", " << g.qubits[1] << ")\n";
        else if (g.name == "swap")
            f << "qc.swap(" << g.qubits[0] << ", " << g.qubits[1] << ")\n";
        else if (g.name == "ccx")
            f << "qc.ccx("  << g.qubits[0] << ", " << g.qubits[1] << ", " << g.qubits[2] << ")\n";
    }

    f << "\nqc.measure_all()\n";
    f << "print(qc.draw())\n";
    std::cout << "Exported Qiskit circuit: " << filename << "\n";
}

void QuantumCircuit::exportCirq(const std::string& filename) const {
    std::ofstream f(filename);
    if (!f.is_open())
        throw std::runtime_error("Cannot open file: " + filename);

    f << "import cirq\n\n";
    f << "q = [cirq.LineQubit(i) for i in range(" << numQubits << ")]\n";
    f << "circuit = cirq.Circuit()\n\n";

    for (const auto& g : circuitLog) {
        if (g.name == "h")
            f << "circuit.append(cirq.H(q[" << g.qubits[0] << "]))\n";
        else if (g.name == "x")
            f << "circuit.append(cirq.X(q[" << g.qubits[0] << "]))\n";
        else if (g.name == "z")
            f << "circuit.append(cirq.Z(q[" << g.qubits[0] << "]))\n";
        else if (g.name == "rx")
            f << "circuit.append(cirq.rx(rads=" << g.param << ")(q[" << g.qubits[0] << "]))\n";
        else if (g.name == "ry")
            f << "circuit.append(cirq.ry(rads=" << g.param << ")(q[" << g.qubits[0] << "]))\n";
        else if (g.name == "rz")
            f << "circuit.append(cirq.rz(rads=" << g.param << ")(q[" << g.qubits[0] << "]))\n";
        else if (g.name == "cx")
            f << "circuit.append(cirq.CNOT(q[" << g.qubits[0] << "], q[" << g.qubits[1] << "]))\n";
        else if (g.name == "swap")
            f << "circuit.append(cirq.SWAP(q[" << g.qubits[0] << "], q[" << g.qubits[1] << "]))\n";
        else if (g.name == "ccx")
            f << "circuit.append(cirq.CCNOT(q[" << g.qubits[0] << "], q[" << g.qubits[1]
              << "], q[" << g.qubits[2] << "]))\n";
    }

    f << "\nprint(circuit)\n";
    f << "simulator = cirq.Simulator()\n";
    f << "result = simulator.simulate(circuit)\n";
    f << "print(result)\n";
    std::cout << "Exported Cirq circuit: " << filename << "\n";
}

void QuantumCircuit::validNoise(double errorRate) {
    if (numQubits > DM_MAX_QUBITS)
        throw std::runtime_error(
            "Noise simulation limited to " + std::to_string(DM_MAX_QUBITS) +
            " qubits (density matrix would require >" + std::to_string(1 << (2 * DM_MAX_QUBITS)) + " complex numbers)");
    noiseModel.depolarizeRate = errorRate;
    noiseModel.active = true;
    dm.initialize(numQubits);
    noiseActivated = true;
}
