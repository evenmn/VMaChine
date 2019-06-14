#pragma once
#include "basis.h"

class HydrogenOrbital : public Basis {
public:
    HydrogenOrbital(System* system);
    //void numberOfOrbitals();
    //void generateListOfStates();

    double evaluate(double x, int n);
    double evaluateDerivative(double x, int n);

    double basisElement(const int n, Eigen::VectorXd positions);
    double basisElementDer(const int n, const int i, Eigen::VectorXd positions);
    double basisElementSecDer(const int n, const int i, Eigen::VectorXd positions);

private:
    int m_Z = 1;
};
