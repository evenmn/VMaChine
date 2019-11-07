#pragma once
#include "hamiltonian.h"

class HarmonicOscillator : public Hamiltonian
{
public:
    HarmonicOscillator(System *system);
    arma::uword getGlobalArrayNeed() { return m_globalArrayNeed; }
    std::string getLabel() { return m_label; }

    void initialize();
    double computeLocalEnergy();
    double getExternalEnergy();

private:
    double m_omega = 0;
    double m_omegaSqrd = 0;
    arma::uword m_globalArrayNeed = 1;
    std::string m_label = "quantumdot";
};
