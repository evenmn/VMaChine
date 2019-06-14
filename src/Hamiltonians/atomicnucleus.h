#pragma once
#include "hamiltonian.h"

class AtomicNucleus : public Hamiltonian {
public:
    AtomicNucleus(System* system);
    double computeLocalEnergy();
    double getExternalEnergy();
    int    getGlobalArrayNeed()  { return m_globalArrayNeed; }
    int    getNumberOfSources()  { return m_numberOfSources; }

private:
    int m_Z = 1;
    int m_globalArrayNeed = 3;
    int m_numberOfSources = 1;
};
