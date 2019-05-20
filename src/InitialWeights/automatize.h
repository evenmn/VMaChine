#pragma once
#include "initialweights.h"

class Automatize : public InitialWeights {
public:
    Automatize(System* system);
    void setupInitialWeights();

private:
    double      m_factor = 0;
    std::string m_trialWaveFunction;
};
