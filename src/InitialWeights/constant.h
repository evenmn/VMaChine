#pragma once
#include "initialweights.h"

class Constant : public InitialWeights
{
public:
    Constant(System *system, const double factor);
    void setupInitialWeights();

    Eigen::MatrixXd getParameters();

private:
    double m_factor = 1;
};
