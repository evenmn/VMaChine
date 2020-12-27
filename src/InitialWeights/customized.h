#pragma once
#include "initialweights.h"

class Customized : public InitialWeights
{
public:
    Customized(System *system);
    void setupInitialWeights();

    Eigen::MatrixXd getParameters();

private:
    double m_factor = 1;
};
