#pragma once
#include "hamiltonian.h"

class AtomicNucleus : public Hamiltonian
{
public:
    AtomicNucleus(System *system);
    double computeLocalEnergy();
    double getExternalEnergy();
    int getGlobalArrayNeed() { return m_globalArrayNeed; }
    std::string getLabel() { return m_label; }

private:
    int m_Z = 1;
    int m_globalArrayNeed = 3;
    std::string m_label = "atom";
};
