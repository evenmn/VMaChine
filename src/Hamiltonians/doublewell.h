#pragma once
#include "hamiltonian.h"

class DoubleWell : public Hamiltonian {
public:
    DoubleWell(System* system, double b);
    double computeLocalEnergy();
    double getExternalEnergy();
    int    getGlobalArrayNeed()  { return m_globalArrayNeed; }
    int    getNumberOfSources()  { return m_numberOfSources; }

private:
    double m_omega = 0;
    double m_omega_sqrd = 0;
    int    m_globalArrayNeed = 1;
    int    m_numberOfSources = 2;
    double m_b = 2;
};
