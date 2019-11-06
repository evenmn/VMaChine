#pragma once
#include "initialweights.h"

class Randomize : public InitialWeights
{
public:
    Randomize(System *system, const double factor = 1);
    void setupInitialWeights();

    arma::mat getParameters();

private:
    double m_factor = 0;
    double m_omega = 1;
};
