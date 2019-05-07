#pragma once
#include "basis.h"

class HydrogenOrbital : public Basis {
public:
    HydrogenOrbital(System* system);
    void            numberOfOrbitals        ();
    void            generateListOfStates    (const int orbitals);

    double          basisElement            (const unsigned int n, const Eigen::VectorXd positions);
    double          basisElementDer         (const unsigned int n, const unsigned int i, const Eigen::VectorXd positions);
    double          basisElementSecDer      (const unsigned int n, const unsigned int i, const Eigen::VectorXd positions);

    double          evaluate                (const double x, const unsigned int n);
    double          evaluateDerivative      (const double x, const unsigned int n);
    double          evaluateSecondDerivative(const double x, const unsigned int n);

private:
    unsigned int    m_Z = 1;
};
