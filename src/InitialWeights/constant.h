#pragma once
#include "initialweights.h"

class Constant : public InitialWeights {
public:
    Constant(System* system, const double factor);
    void setupInitialWeights();

private:
    double m_factor = 1;
};
