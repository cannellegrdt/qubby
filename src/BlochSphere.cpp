/*
 * Project: qubby
 * File name: BlochSphere.cpp
 * Author: Cannelle Gourdet - lankley
 * File description: Implements BlochSphere: Bloch vector computation, ASCII terminal art,
 *                   and SVG export for single-qubit state visualisation.
 */

#include "BlochSphere.hpp"
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <algorithm>

/**
 * @details
 * Computes the reduced density matrix ρ_q for qubit @p qubit by summing over all
 * other qubit degrees of freedom. For each pair (s, s|bit) where bit = (1<<qubit):
 *
 *   ρ_{00} += |ψ_s|²         (s has qubit q = |0⟩)
 *   ρ_{11} += |ψ_{s|bit}|²   (index with qubit q flipped to |1⟩)
 *   ρ_{01} += ψ_s · conj(ψ_{s|bit})
 *
 * Bloch vector:
 *   x = Tr(X ρ_q) = 2 Re(ρ_{01})
 *   y = Tr(Y ρ_q) = −2 Im(ρ_{01})
 *   z = Tr(Z ρ_q) = ρ_{00} − ρ_{11}
 */
BlochSphere::BlochVector BlochSphere::compute(const QuantumState& state, int qubit) {
    int n = state.getNumQubits();
    if (qubit < 0 || qubit >= n)
        throw std::runtime_error("Bloch sphere: qubit index " + std::to_string(qubit) +
                                 " out of range [0, " + std::to_string(n - 1) + "]");

    int dim = 1 << n;
    int bit = 1 << qubit;

    std::complex<double> rho01 = 0.0;
    double z = 0.0;

    for (int i = 0; i < dim; i++) {
        if ((i & bit) == 0) {
            auto a = state.getAmplitude(i);
            auto b = state.getAmplitude(i | bit);
            rho01 += a * std::conj(b);
            z += std::norm(a) - std::norm(b);
        }
    }

    double bx = 2.0 * rho01.real();
    double by = -2.0 * rho01.imag();
    double bz = z;
    double norm = std::sqrt(bx * bx + by * by + bz * bz);

    return {bx, by, bz, norm};
}

/**
 * @details
 * Renders the sphere as a 35-column × 19-row character grid (xz-plane projection).
 * Character aspect ratio compensation: horizontal half = 15 chars, vertical half = 8 chars.
 *
 * Draw order (later writes overwrite earlier):
 *   1. Circle outline ('o') via parametric scan (360 steps).
 *   2. z-axis ('|') and x-axis ('-'); intersection ('+').
 *   3. State vector line ('.') from centre using linear interpolation.
 *   4. State vector endpoint ('*').
 */
void BlochSphere::printASCII(const BlochVector& bv, int qubit) {
    static const int W = 35;
    static const int H = 19;
    static const int cx = 17;
    static const int cz = 9;
    static const double hw = 14.5;
    static const double hh =  8.0;

    std::vector<std::string> grid(H, std::string(W, ' '));

    for (int deg = 0; deg < 360; deg++) {
        double t = deg * M_PI / 180.0;
        int col = cx + (int)std::round(hw * std::cos(t));
        int row = cz - (int)std::round(hh * std::sin(t));
        if (row >= 0 && row < H && col >= 0 && col < W)
            grid[row][col] = 'o';
    }

    for (int r = 0; r < H; r++) grid[r][cx] = '|';
    for (int c = 0; c < W; c++) grid[cz][c]  = '-';
    grid[cz][cx] = '+';

    int vcol = cx + (int)std::round(hw * bv.x);
    int vrow = cz - (int)std::round(hh * bv.z);
    vcol = std::clamp(vcol, 0, W - 1);
    vrow = std::clamp(vrow, 0, H - 1);

    if (bv.norm > 0.02) {
        int dc = vcol - cx, dr = vrow - cz;
        int steps = std::max(std::abs(dc), std::abs(dr));
        for (int t = 1; t < steps; t++) {
            int c = cx + dc * t / steps;
            int r = cz + dr * t / steps;
            if (grid[r][c] == ' ')
                grid[r][c] = '.';
        }
        grid[vrow][vcol] = '*';
    }

    std::string title = "  Bloch sphere — q[" + std::to_string(qubit) + "]  (xz projection)";
    std::cout << "\n" << title << "\n\n";
    std::cout << "       |0" << "\xe2\x9f\xa9" << "\n";
    for (int r = 0; r < H; r++) {
        if (r == cz)
            std::cout << "-x " << grid[r] << " +x\n";
        else
            std::cout << "   " << grid[r] << "\n";
    }
    std::cout << "       |1" << "\xe2\x9f\xa9" << "\n\n";

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "  r = (" << bv.x << ", " << bv.y << ", " << bv.z << ")"
              << "   |r| = " << bv.norm;
    if (bv.norm > 0.999)
        std::cout << " (pure)";
    else if (bv.norm < 0.01)
        std::cout << " (fully mixed)";
    else
        std::cout << " (mixed — entangled)";
    std::cout << "\n  y-component not shown in terminal view.\n\n";
}

namespace {

/* Orthographic-ish projection used throughout the SVG.
 * y-axis goes toward upper-right (ky, kzy controls).
 * z-axis goes up, x-axis goes right.
 *
 *  sx = cx_svg + r_svg * (bx + ky  * by)
 *  sy = cy_svg - r_svg * (bz + kzy * by)
 */
struct Proj {
    double cx, cy, r;
    static constexpr double ky = 0.28;
    static constexpr double kzy = 0.17;

    double sx(double bx, double by) const { return cx + r * (bx + ky  * by); }
    double sy(double bz, double by) const { return cy - r * (bz + kzy * by); }
};

std::string equatorPath(const Proj& p) {
    std::ostringstream ss;
    for (int i = 0; i <= 64; i++) {
        double t  = 2.0 * M_PI * i / 64.0;
        double ex = p.sx(std::cos(t), std::sin(t));
        double ey = p.sy(0.0, std::sin(t));
        ss << std::fixed << std::setprecision(1);
        if (i == 0) ss << "M " << ex << " " << ey;
        else ss << " L " << ex << " " << ey;
    }
    return ss.str();
}

/* Axis endpoint projected to sphere surface (unit vector). */
struct AxisPt {
    double x, y;
};
AxisPt axisEnd(const Proj& p, double bx, double by, double bz) {
    return {p.sx(bx, by), p.sy(bz, by)};
}

} // anonymous namespace

/**
 * @details
 * SVG layout (420 × 450):
 *   - Centre of sphere at (210, 215), radius 160.
 *   - Projection: orthographic with y-axis tilted upper-right (28% horizontal, 17% vertical).
 *   - Defs: two SVG markers for axis and vector arrowheads.
 *   - Elements: sphere circle, dashed equatorial ellipse, negative-axis dashed lines,
 *     positive-axis lines with arrowheads, state vector arrow, labels, coordinates.
 *
 * @note Uses explicit stream writes (no raw-string literals) to avoid the conflict
 *       between the `)"` raw-string terminator and `url(#...)` SVG references.
 */
void BlochSphere::saveSVG(const BlochVector& bv, int qubit, const std::string& filename) {
    std::ofstream f(filename);
    if (!f.is_open())
        throw std::runtime_error("Cannot open SVG file: " + filename);

    Proj p{210.0, 215.0, 160.0};

    auto pz_pos = axisEnd(p, 0, 0, 1);
    auto pz_neg = axisEnd(p, 0, 0, -1);
    auto px_pos = axisEnd(p, 1, 0, 0);
    auto px_neg = axisEnd(p, -1, 0, 0);
    auto py_pos = axisEnd(p, 0, 1, 0);
    auto py_neg = axisEnd(p, 0, -1, 0);

    double vsx = p.sx(bv.x, bv.y);
    double vsy = p.sy(bv.z, bv.y);

    auto d1 = [](double v) {
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(1) << v;
        return ss.str();
    };
    auto d4 = [](double v) {
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(4) << v;
        return ss.str();
    };
    auto line = [&](double x1, double y1, double x2, double y2,
                    const std::string& cls, const std::string& extra = "") {
        f << "<line x1=\"" << d1(x1) << "\" y1=\"" << d1(y1)
          << "\" x2=\"" << d1(x2) << "\" y2=\"" << d1(y2)
          << "\" class=\"" << cls << "\"" << extra << "/>\n";
    };

    f << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
      << "<svg width=\"420\" height=\"450\" viewBox=\"0 0 420 450\""
      << " xmlns=\"http://www.w3.org/2000/svg\">\n"
      << "<defs>\n"
      << "  <marker id=\"axis-arr\" viewBox=\"0 0 10 10\" refX=\"9\" refY=\"5\"\n"
      << "          markerWidth=\"6\" markerHeight=\"6\" orient=\"auto\">\n"
      << "    <path d=\"M 0 0 L 10 5 L 0 10 z\" fill=\"#444\"/>\n"
      << "  </marker>\n"
      << "  <marker id=\"vec-arr\" viewBox=\"0 0 10 10\" refX=\"9\" refY=\"5\"\n"
      << "          markerWidth=\"7\" markerHeight=\"7\" orient=\"auto\">\n"
      << "    <path d=\"M 0 0 L 10 5 L 0 10 z\" fill=\"#cc2222\"/>\n"
      << "  </marker>\n"
      << "</defs>\n"
      << "<style>\n"
      << "  .lbl   { font: 14px sans-serif; fill: #333; }\n"
      << "  .coord { font: 12px monospace;  fill: #555; }\n"
      << "  .axis  { stroke: #444; stroke-width: 1.5; fill: none; }\n"
      << "  .dash  { stroke: #aaa; stroke-width: 1;   fill: none; stroke-dasharray: 5,4; }\n"
      << "  .equator { stroke: #888; stroke-width: 1; fill: none; stroke-dasharray: 6,4; }\n"
      << "  .sphere  { stroke: #333; stroke-width: 1.5; fill: #f4f4fb; }\n"
      << "  .vector  { stroke: #cc2222; stroke-width: 2.5; fill: none; }\n"
      << "  .vdot    { fill: #cc2222; }\n"
      << "</style>\n\n";

    f << "<!-- sphere -->\n"
      << "<circle cx=\"" << d1(p.cx) << "\" cy=\"" << d1(p.cy)
      << "\" r=\"" << d1(p.r) << "\" class=\"sphere\"/>\n\n";

    f << "<!-- equatorial ellipse -->\n"
      << "<path d=\"" << equatorPath(p) << "\" class=\"equator\"/>\n\n";

    f << "<!-- negative axes -->\n";
    line(p.cx, p.cy, pz_neg.x, pz_neg.y, "dash");
    line(p.cx, p.cy, px_neg.x, px_neg.y, "dash");
    line(p.cx, p.cy, py_neg.x, py_neg.y, "dash");

    const std::string aa = " marker-end=\"url(#axis-arr)\"";
    f << "\n<!-- positive axes -->\n";
    line(p.cx, p.cy, pz_pos.x, pz_pos.y, "axis", aa);
    line(p.cx, p.cy, px_pos.x, px_pos.y, "axis", aa);
    line(p.cx, p.cy, py_pos.x, py_pos.y, "axis", aa);

    f << "\n<!-- axis labels -->\n"
      << "<text x=\"" << d1(pz_pos.x) << "\" y=\"" << d1(pz_pos.y - 12)
      << "\" text-anchor=\"middle\" class=\"lbl\">|0&#x27E9;</text>\n"
      << "<text x=\"" << d1(pz_neg.x) << "\" y=\"" << d1(pz_neg.y + 20)
      << "\" text-anchor=\"middle\" class=\"lbl\">|1&#x27E9;</text>\n"
      << "<text x=\"" << d1(px_pos.x + 10) << "\" y=\"" << d1(px_pos.y + 5)
      << "\" class=\"lbl\">x</text>\n"
      << "<text x=\"" << d1(py_pos.x + 6) << "\" y=\"" << d1(py_pos.y - 6)
      << "\" class=\"lbl\">y</text>\n\n";

    const std::string va = " marker-end=\"url(#vec-arr)\"";
    if (bv.norm > 0.01) {
        f << "<!-- state vector -->\n";
        line(p.cx, p.cy, vsx, vsy, "vector", va);
        f << "<circle cx=\"" << d1(vsx) << "\" cy=\"" << d1(vsy)
          << "\" r=\"5\" class=\"vdot\"/>\n\n";
    } else {
        f << "<!-- state vector at origin (fully mixed) -->\n"
          << "<circle cx=\"" << d1(p.cx) << "\" cy=\"" << d1(p.cy)
          << "\" r=\"6\" class=\"vdot\"/>\n\n";
    }

    std::string purity = (bv.norm > 0.999) ? "pure" :
                         (bv.norm < 0.01) ? "fully mixed" : "mixed";
    f << "<!-- coordinates -->\n"
      << "<text x=\"10\" y=\"430\" class=\"coord\">q[" << qubit << "]  "
      << "x=" << d4(bv.x) << "  y=" << d4(bv.y) << "  z=" << d4(bv.z)
      << "  |r|=" << d4(bv.norm) << "  (" << purity << ")</text>\n\n"
      << "</svg>\n";

    std::cout << "Saved Bloch sphere SVG: " << filename << "\n";
}
