#pragma once
#include "basis.h"

class Hermite : public Basis {
public:
    Hermite();
    double evaluate(double x, int n);
    double evaluateDerivative(double x, int n);
};
