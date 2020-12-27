#pragma once
#include "hamiltonian.h"

class DoubleWell : public Hamiltonian
{
public:
    DoubleWell(System *system, double displacement = 2);
    int getGlobalArrayNeed() { return m_globalArrayNeed; }
    std::string getLabel() { return m_label; }

    void initialize();
    double computeLocalEnergy();
    double getExternalEnergy();

private:
    double m_omega = 0;
    double m_omega_sqrd = 0;
    int m_globalArrayNeed = 1;
    double m_b = 2;
    double m_offset = 1;
    std::string m_label = "doubledot";
};
