#pragma once
#include "hamiltonian.h"
#include "../Eigen/Dense"

class AtomicNucleus : public Hamiltonian {
public:
    AtomicNucleus(System* system);
    double computeLocalEnergy();
    double getExternalEnergy();
    int    getGlobalArrayNeed()  { return m_globalArrayNeed; }

private:
    int m_Z = 1;
    int m_globalArrayNeed = 3;
};
