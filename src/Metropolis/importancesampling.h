#pragma once
#include "metropolis.h"
#include <Eigen/Dense>
#include <vector>

class ImportanceSampling : public Metropolis {
public:
    ImportanceSampling(System* system);
    bool acceptMove();
    double QuantumForce(const int i);
    double GreenFuncSum();

    double calculateDistanceMatrixElement(int i, int j);
    void   calculateDistanceMatrixCross(int particle);
    double calculateRadialVectorElement(int particle);

private:
    Eigen::VectorXd m_quantumForceOld;
    Eigen::VectorXd m_quantumForceNew;
    std::vector<class WaveFunction*>    m_waveFunctionVector;
};
