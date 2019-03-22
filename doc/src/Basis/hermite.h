#pragma once
#include "basis.h"

class Hermite : public Basis {
public:
    Hermite(System* system);
    void numberOfOrbitals();
    double evaluate(double x, int n);
    double evaluateDerivative(double x, int n);
    double evaluateSecondDerivative(double x, int n);
    Eigen::MatrixXd generateListOfStates();

private:
    double m_omega = 1;
};
