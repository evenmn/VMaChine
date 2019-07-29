#pragma once
#include "metropolis.h"
#include <vector>

class ImportanceSampling : public Metropolis
{
public:
    ImportanceSampling(System *system);
    bool acceptMove();
    double QuantumForce(const int i);
    double GreenFuncSum();

private:
    Eigen::VectorXd m_quantumForceOld;
    Eigen::VectorXd m_quantumForceNew;
    std::vector<class WaveFunction *> m_waveFunctionVector;

    double m_dtD = 1;
    double m_sqrtStep = 0;
};
