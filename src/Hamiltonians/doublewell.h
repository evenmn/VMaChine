#pragma once
#include "hamiltonian.h"
#include <Eigen/Dense>

class DoubleWell : public Hamiltonian {
public:
    DoubleWell(System* system);
    double          computeLocalEnergy();
    unsigned int    getGlobalArrayNeed()  { return m_globalArrayNeed; }

private:
    unsigned int    m_globalArrayNeed   = 1;
    signed   int    m_A                 = 2;

    double          m_omega             = 0;
    double          m_omega_sqrd        = 0;
};
