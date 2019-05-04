#pragma once
#include "hamiltonian.h"
#include <Eigen/Dense>

class DoubleWell : public Hamiltonian {
public:
    DoubleWell(System* system);
    double computeLocalEnergy();
    int    getGlobalArrayNeed()  { return m_globalArrayNeed; }

private:
    double m_omega = 0;
    double m_omega_sqrd = 0;
    int    m_globalArrayNeed = 1;
};
