#pragma once
#include <Eigen/Dense>

class Basis {
public:
    Basis();
    virtual double evaluate(double x, int n) = 0;
    virtual double evaluateDerivative(double x, int n) = 0;
    virtual ~Basis() = 0;
};
