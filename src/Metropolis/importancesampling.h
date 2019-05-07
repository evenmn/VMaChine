#pragma once
#include "metropolis.h"
#include <Eigen/Dense>
#include <vector>

class ImportanceSampling : public Metropolis {
public:
    ImportanceSampling(System* system);
    bool acceptMove();

    double QuantumForce(const unsigned int i);
    double GreenFuncSum();

private:
    Eigen::VectorXd m_quantumForceOld;
    Eigen::VectorXd m_quantumForceNew;

    std::vector<class WaveFunction*>    m_waveFunctionVector;
};
