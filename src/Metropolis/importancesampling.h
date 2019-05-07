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

    double calculateDistanceMatrixElement   (const unsigned int i, const unsigned int j);
    void   calculateDistanceMatrixCross     (const unsigned int particle);
    double calculateRadialVectorElement     (const unsigned int particle);

private:
    Eigen::VectorXd m_quantumForceOld;
    Eigen::VectorXd m_quantumForceNew;

    std::vector<class WaveFunction*>    m_waveFunctionVector;
};
