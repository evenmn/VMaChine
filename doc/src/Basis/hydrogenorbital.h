#pragma once
#include "basis.h"

class HydrogenOrbital : public Basis {
public:
    HydrogenOrbital(System* system);
    void numberOfOrbitals();
    double evaluate(double x, int n);
    double evaluateDerivative(double x, int n);
    Eigen::MatrixXd generateListOfStates();

private:
    int m_Z = 1;
};
