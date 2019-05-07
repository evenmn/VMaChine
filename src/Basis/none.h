#pragma once
#include "basis.h"

class None : public Basis {
public:
    None(System* system);
    void    numberOfOrbitals();
    void    generateListOfStates(const int orbitals);

    double  basisElement        (const unsigned int n, const Eigen::VectorXd positions);
    double  basisElementDer     (const unsigned int n, const unsigned int i, const Eigen::VectorXd positions);
    double  basisElementSecDer  (const unsigned int n, const unsigned int i, const Eigen::VectorXd positions);
};
