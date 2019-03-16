#pragma once
#include "basis.h"

class None : public Basis {
public:
    None(System* system);
    void numberOfOrbitals();
    double evaluate(double x, int n);
    double evaluateDerivative(double x, int n);
    Eigen::MatrixXd generateListOfStates();
};
