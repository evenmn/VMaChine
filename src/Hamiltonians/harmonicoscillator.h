#pragma once
#include "hamiltonian.h"
#include <Eigen/Dense>

class HarmonicOscillator : public Hamiltonian {
public:
    HarmonicOscillator(System* system);
    double          computeLocalEnergy();
    unsigned int    getGlobalArrayNeed()  { return m_globalArrayNeed; }

private:
    unsigned int    m_globalArrayNeed = 1;
    double          m_omega           = 0;
    double          m_omega_sqrd      = 0;
};

