#pragma once
#include "hamiltonian.h"

class HarmonicOscillator : public Hamiltonian {
public:
    HarmonicOscillator(System* system);
    double computeLocalEnergy();
    double getExternalEnergy();
    int    getGlobalArrayNeed()  { return m_globalArrayNeed; }
    int    getNumberOfSources()  { return m_numberOfSources; }

private:
    double m_omega = 0;
    double m_omegaSqrd = 0;
    int    m_globalArrayNeed = 1;
    int    m_numberOfSources = 1;
};

