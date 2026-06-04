#ifndef __qiskitcpp_primitives_backend_estimator_v2_hpp__
#define __qiskitcpp_primitives_backend_estimator_v2_hpp__

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "circuit/quantumcircuit.hpp"
#include "compiler/transpiler.hpp"
#include "primitives/backend_sampler_v2.hpp"
#include "primitives/containers/sampler_pub.hpp"
#include "providers/backend.hpp"

namespace Qiskit {
namespace primitives {

struct EstimatorPub {
    circuit::QuantumCircuit circuit;
    std::vector<std::pair<std::string, double>> observables;
    uint_t shots = 0;

    EstimatorPub() {}

    EstimatorPub(
        circuit::QuantumCircuit& circ,
        std::vector<std::pair<std::string, double>> obs,
        uint_t in_shots = 0
    ) : circuit(circ), observables(obs), shots(in_shots) {}
};

struct EstimatorPubResult {
    double ev = 0.0;
    double stddev = 0.0;
};

class BackendEstimatorV2 {
protected:
    uint_t shots_;
    providers::BackendV2& backend_;

    static int bit_value_from_string(
        const std::string& bitstring,
        uint_t qubit,
        uint_t num_qubits
    ) {
        uint_t pos = num_qubits - 1 - qubit;
        return bitstring[pos] == '1' ? 1 : 0;
    }

    static double pauli_eigenvalue_from_bitstring(
        const std::string& pauli,
        const std::string& bitstring
    ) {
        double value = 1.0;
        uint_t n = pauli.size();

        for (uint_t q = 0; q < n; q++) {
            char p = pauli[n - 1 - q];

            if (p == 'I') {
                continue;
            }

            int bit = bit_value_from_string(bitstring, q, n);
            value *= (bit == 0) ? 1.0 : -1.0;
        }

        return value;
    }

    static void append_basis_rotation(
        circuit::QuantumCircuit& circ,
        const std::string& pauli
    ) {
        uint_t n = pauli.size();

        for (uint_t q = 0; q < n; q++) {
            char p = pauli[n - 1 - q];

            if (p == 'X') {
                circ.h(q);
            } else if (p == 'Y') {
                circ.rz(-M_PI / 2.0, q);
                circ.h(q);
            }
        }
    }

    static void append_measurements(circuit::QuantumCircuit& circ, uint_t num_qubits) {
        for (uint_t q = 0; q < num_qubits; q++) {
            circ.measure(q, q);
        }
    }

public:
    BackendEstimatorV2(providers::BackendV2& backend, uint_t shots = 1024)
        : shots_(shots), backend_(backend) {}

    const providers::BackendV2& backend(void) const {
        return backend_;
    }

    std::vector<EstimatorPubResult> run(std::vector<EstimatorPub> pubs) {
        std::vector<EstimatorPubResult> results;
        BackendSamplerV2 sampler(backend_, shots_);

        for (auto& pub : pubs) {
            double total_ev = 0.0;
            double variance_accum = 0.0;

            for (auto& term : pub.observables) {
                const std::string& pauli = term.first;
                double coeff = term.second;
                uint_t num_qubits = pauli.size();

                auto meas_circ = pub.circuit.copy();

                append_basis_rotation(meas_circ, pauli);
                append_measurements(meas_circ, num_qubits);

                auto isa_circ = compiler::transpile(meas_circ, backend_);

                SamplerPub sampler_pub(
                    isa_circ,
                    pub.shots > 0 ? pub.shots : shots_
                );

                auto job = sampler.run({sampler_pub});
                auto primitive_result = job->result();

                auto pub_result = primitive_result[0];

                auto bits = pub_result.data().get_bitstrings();

                if (bits.empty()) {
                    continue;
                }

                double mean = 0.0;
                for (auto& b : bits) {
                    mean += pauli_eigenvalue_from_bitstring(pauli, b);
                }

                mean /= static_cast<double>(bits.size());

                total_ev += coeff * mean;

                double var =
                    (1.0 - mean * mean) / static_cast<double>(bits.size());

                variance_accum += coeff * coeff * var;
            }

            EstimatorPubResult r;
            r.ev = total_ev;
            r.stddev = std::sqrt(variance_accum);
            results.push_back(r);
        }

        return results;
    }
};

} // namespace primitives
} // namespace Qiskit

#endif
