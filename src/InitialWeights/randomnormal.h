#pragma once
#include "initialweights.h"

class RandomNormalWeights : public InitialWeights
{
public:
    RandomNormalWeights(System *system, const double variance = 1);
    void setupInitialWeights();

    Eigen::MatrixXd getParameters();

private:
    double m_variance = 0;
};
