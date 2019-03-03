#pragma once
#include "basis.h"

class HydrogenOrbital : public Basis {
public:
    HydrogenOrbital(System* system, int Z);
    void numberOfOrbitals();
    double evaluate(double x, int n);
    double evaluateDerivative(double x, int n);

private:
    int m_Z = 1;
};
