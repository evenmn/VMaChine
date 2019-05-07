#pragma once
#include "hamiltonian.h"
#include <Eigen/Dense>

class AtomicNucleus : public Hamiltonian {
public:
    AtomicNucleus(System* system);
    double          computeLocalEnergy();
    unsigned int    getGlobalArrayNeed()  { return m_globalArrayNeed; }

private:
    unsigned int m_Z                = 1;
    unsigned int m_globalArrayNeed  = 3;
};
