/*
 * Project: qubby
 * File name: NoiseModel.hpp
 * Author: Cannelle Gourdet - lankley
 * File description: Defines the NoiseModel struct, which holds per-gate error rates
 *                   for the three standard quantum noise channels.
 */

#ifndef NOISEMODEL_HPP_
    #define NOISEMODEL_HPP_

/**
 * @struct NoiseModel
 * @brief Parameters for realistic noise simulation via quantum channels.
 *
 * When active, the noise model is applied after every gate in QuantumCircuit::load().
 * Each rate is the probability that the corresponding error occurs on a qubit per gate:
 *
 *   - bitFlipRate    → bit-flip   (X error)   channel
 *   - phaseFlipRate  → phase-flip (Z error)   channel
 *   - depolarizeRate → depolarizing (X+Y+Z)  channel
 *
 * Rates must be in [0, 1]. Multiple channels stack independently.
 */
struct NoiseModel {
    double bitFlipRate = 0.0;
    double phaseFlipRate = 0.0;
    double depolarizeRate = 0.0;
    bool active = false;
};

#endif /* NOISEMODEL_HPP_ */
