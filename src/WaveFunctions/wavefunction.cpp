#include "wavefunction.h"
#include <iostream>

WaveFunction::WaveFunction(System * system) {
    m_system = system;
}

Eigen::Map<Eigen::VectorXd> WaveFunction::flatten(Eigen::MatrixXd A) {
    return Eigen::Map<Eigen::VectorXd>(A.data(), A.size());
}

Eigen::Map<Eigen::MatrixXd> WaveFunction::reshape(Eigen::VectorXd A, const Eigen::Index m, const Eigen::Index n) {
    return Eigen::Map<Eigen::MatrixXd>(A.data(), m, n);
}

WaveFunction::~WaveFunction() {}
