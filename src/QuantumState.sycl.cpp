/*
 * Project: qubby
 * File name: QuantumState.sycl.cpp
 * Author: Cannelle Gourdet - lankley
 * File description: SYCL implementation of QuantumState::applyGateSYCL(), offloading
 *                   the gate application loop to the GPU via Intel oneAPI.
 */

#include "QuantumState.hpp"
#include <sycl/sycl.hpp>
#include <array>

#ifdef SYCL_ENABLED
    /**
     * @details
     * Data flow:
     * 1. The 2×2 gate is flattened into a 4-element `std::array` (row-major: m00, m01, m10, m11)
     *    so it can be transferred to the device as a contiguous SYCL buffer.
     * 2. Two SYCL buffers are created — one for the full amplitude vector (read/write),
     *    one for the gate coefficients (read-only). Buffer destruction triggers implicit
     *    host←device copy-back, so no explicit memcpy is needed.
     * 3. A single `parallel_for` kernel is submitted: each work-item idx handles
     *    one half-index i = idx[0], reconstructs i0 and i1 via the same bitwise
     *    insertion as applyGate(), reads val0/val1 before writing to avoid races,
     *    then writes the updated amplitudes.
     * 4. `q.wait()` blocks until all work-items have finished before returning.
     *
     * @note Race-free by construction: each (i0, i1) pair is handled by exactly
     *       one work-item. No barriers or atomics are needed.
     */
    void QuantumState::applyGateSYCL(int targetQubit, Matrix quantumGate) {
        size_t loop_size = 1ULL << (this->num_qubits - 1);
        std::array<std::complex<double>, 4> gate_data = {
            quantumGate.m00, quantumGate.m01,
            quantumGate.m10, quantumGate.m11
        };

        sycl::queue q(sycl::default_selector_v);

        sycl::buffer<std::complex<double>, 1> buffer_amps(this->amplitudes.data(), sycl::range<1>(amplitudes.size()));
        sycl::buffer<std::complex<double>, 1> buffer_gate(gate_data.data(), sycl::range<1>(4));

        q.submit([&](sycl::handler& h) {
            auto acc_amps = buffer_amps.get_access<sycl::access::mode::read_write>(h);
            auto acc_gate = buffer_gate.get_access<sycl::access::mode::read>(h);

            h.parallel_for(sycl::range<1>(loop_size), [=](sycl::id<1> idx) {
                size_t i = idx[0];

                uint64_t mask = (1ULL << targetQubit) - 1;
                uint64_t i0 = ((i >> targetQubit) << (targetQubit + 1)) | (i & mask);
                uint64_t i1 = i0 | (1ULL << targetQubit);

                std::complex<double> val0 = acc_amps[i0];
                std::complex<double> val1 = acc_amps[i1];

                acc_amps[i0] = acc_gate[0] * val0 + acc_gate[1] * val1;
                acc_amps[i1] = acc_gate[2] * val0 + acc_gate[3] * val1;
            });
        });

        q.wait();
    }
#endif
