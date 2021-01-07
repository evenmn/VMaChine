#pragma once
#include "initialweights.h"

class RandomUniformWeights : public InitialWeights
{
public:
    RandomUniformWeights(System *system, const double factor = 1);
    std::string getLabel() { return m_label; }
    void setupInitialWeights();

    Eigen::MatrixXd getParameters();

private:
    double m_factor = 0;
    std::string m_label = "random uniform";
};
