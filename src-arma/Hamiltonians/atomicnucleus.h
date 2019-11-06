#pragma once
#include "hamiltonian.h"

class AtomicNucleus : public Hamiltonian
{
public:
    AtomicNucleus(System *system);
    arma::uword getGlobalArrayNeed() { return m_globalArrayNeed; }
    std::string getLabel() { return m_label; }

    void initialize();
    double computeLocalEnergy();
    double getExternalEnergy();

private:
    arma::uword m_Z = 1;
    arma::uword m_globalArrayNeed = 3;
    std::string m_label = "atom";
};
