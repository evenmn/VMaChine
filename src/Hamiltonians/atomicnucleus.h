#pragma once
#include "hamiltonian.h"

class AtomicNucleus : public Hamiltonian
{
public:
    AtomicNucleus(System *system);
    int getGlobalArrayNeed() { return m_globalArrayNeed; }
    std::string getLabel() { return m_label; }

    void initialize();
    double getExternalEnergy();

private:
    int m_Z = 1;
    int m_globalArrayNeed = 3;
    std::string m_label = "atom";
};
