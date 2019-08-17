#pragma once
#include "hamiltonian.h"

class HarmonicOscillator : public Hamiltonian
{
public:
    HarmonicOscillator(System *system);
    double computeLocalEnergy();
    double getExternalEnergy();
    int getGlobalArrayNeed() { return m_globalArrayNeed; }
    std::string getLabel() { return m_label; }

private:
    double m_omega = 0;
    double m_omegaSqrd = 0;
    int m_globalArrayNeed = 1;
    std::string m_label = "quantumdot";
};
