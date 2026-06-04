/*
# This code is part of Qiskit.
#
# (C) Copyright IBM 2026.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.
*/

#include <iostream>
#include <vector>
#include <utility>
#include <string>
#include <cstdlib>

#include "circuit/quantumcircuit.hpp"
#include "primitives/backend_estimator_v2.hpp"
#include "service/qiskit_runtime_service_c.hpp"
#include "compiler/transpiler.hpp"

using namespace Qiskit;
using namespace Qiskit::circuit;
using namespace Qiskit::primitives;
using namespace Qiskit::service;
using namespace Qiskit::compiler;

using Estimator = BackendEstimatorV2;

int main() {
    int num_qubits = 2;

    auto qreg = QuantumRegister(num_qubits);
    auto creg = ClassicalRegister(num_qubits, std::string("meas"));
    QuantumCircuit circ(
        std::vector<QuantumRegister>({qreg}),
        std::vector<ClassicalRegister>({creg})
    );

    // Bell state
    circ.h(0);
    circ.cx(0, 1);
    circ.measure(qreg, creg);

    // Use your environment/saved account
    auto service = QiskitRuntimeService();

    // Replace with an available backend/simulator
    auto backend = service.backend("ibm_kingston");

    auto transpiled_circ = transpile(circ, backend);

    std::vector<std::pair<std::string, double>> obs = {
        {"ZZ", 1.0},
        {"XI", 0.5}
    };

    auto estimator = Estimator(backend, 1000);
    auto results = estimator.run({EstimatorPub(transpiled_circ, obs, 1000)});

    std::cout << "EV = " << results[0].ev
              << " +/- " << results[0].stddev << std::endl;

    return 0;
}