#pragma once
#include "hamiltonian.h"

class EllipticalHarmonicOscillator : public Hamiltonian
{
public:
    EllipticalHarmonicOscillator(System *system, const double beta);
    int getGlobalArrayNeed() { return m_globalArrayNeed; }
    std::string getLabel() { return m_label; }

    void initialize();
    double getExternalEnergy();

private:
    double m_omegaSqrd = 0;
    double m_beta = sqrt(8);
    int m_globalArrayNeed = 1;
    std::string m_label = "elliptical harmonic oscillator";
};
