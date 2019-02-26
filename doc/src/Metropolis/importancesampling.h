#pragma once
#include "metropolis.h"
#include <Eigen/Dense>
#include <vector>

class ImportanceSampling : public Metropolis {
public:
    ImportanceSampling(System* system);
    bool acceptMove();
    double QuantumForce(const int i);
    double GreenFuncSum(const Eigen::VectorXd newPositions, const Eigen::VectorXd QuantumForceOld);

private:
    int m_numberOfParticles = 0;
    int m_numberOfDimensions = 0;
    std::vector<class WaveFunction*>    m_waveFunctionVector;
};
