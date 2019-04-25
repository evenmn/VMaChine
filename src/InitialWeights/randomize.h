#pragma once
#include "initialweights.h"

class Randomize : public InitialWeights {
public:
    Randomize(System* system, const double factor);
    void setupInitialWeights();

private:
    double m_factor = 0;
};