#pragma once
#include "initialweights.h"

class RandomUniformWeights : public InitialWeights
{
public:
    RandomUniformWeights(System *system, const double factor = 1);
    void setupInitialWeights();

    Eigen::MatrixXd getParameters();

private:
    double m_factor = 0;
};
