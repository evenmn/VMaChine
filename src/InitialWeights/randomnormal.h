#pragma once
#include "initialweights.h"

class RandomNormalWeights : public InitialWeights
{
public:
    RandomNormalWeights(System *system, const double variance = 1);
    std::string getLabel() { return m_label; }
    void setupInitialWeights();

    Eigen::MatrixXd getParameters();

private:
    double m_variance = 0;
    std::string m_label = "random normal";
};
