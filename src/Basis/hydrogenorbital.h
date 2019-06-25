#pragma once
#include "basis.h"

class HydrogenOrbital : public Basis {
public:
    HydrogenOrbital(System* system);

    double evaluate(double x, int n);
    double evaluateDerivative(double x, int n);
    double evaluateSecondDerivative(double x, int n);

    double basisElement(const int n, Eigen::VectorXd positions);
    double basisElementDer(const int n, const int i, Eigen::VectorXd positions);
    double basisElementSecDer(const int n, const int i, Eigen::VectorXd positions);

    double radial(double r, int n);
    double angular(double theta, double phi, int l, int m);

private:
    int m_Z = 1;
};
