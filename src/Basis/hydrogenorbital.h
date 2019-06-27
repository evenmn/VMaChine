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

    double evaluateCart(Eigen::VectorXd position, int n, int l, int m);
    double evaluateCartDerivative(Eigen::VectorXd position, int i, int n, int l, int m);
    double evaluateCartSecondDerivative(Eigen::VectorXd position, int i, int n, int l, int m);

private:
    int m_Z = 1;
    double m_alpha = 1;
};
