#pragma once
#include "hamiltonian.h"
#include <Eigen/Dense>

class AtomicNucleus : public Hamiltonian {
public:
    AtomicNucleus(System* system);
    double computeLocalEnergy();

private:
    int m_Z = 1;
};
