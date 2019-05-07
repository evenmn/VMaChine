#pragma once
#include "basis.h"

class Hermite : public Basis {
public:
    Hermite(System* system);
    void            numberOfOrbitals         ();
    void            generateListOfStates     (const int orbitals);

    double          basisElement             (const unsigned int n, const Eigen::VectorXd positions);
    double          basisElementDer          (const unsigned int n, const unsigned int i, Eigen::VectorXd positions);
    double          basisElementSecDer       (const unsigned int n, const unsigned int i, Eigen::VectorXd positions);

    double          evaluate                 (const double x, const int n);
    double          evaluateDerivative       (const double x, const int n);
    double          evaluateSecondDerivative (const double x, const int n);

private:
    double          m_omega         = 1;
    double          m_omegaSqrt     = 1;
    Eigen::MatrixXi m_listOfStates;
};
