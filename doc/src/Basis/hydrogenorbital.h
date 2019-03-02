#pragma once
#include "basis.h"

class HydrogenOrbital : public Basis {
public:
    HydrogenOrbital(int Z);
    double evaluate(double x, int n);
    double evaluateDerivative(double x, int n);

private:
    int m_Z = 1;
};
