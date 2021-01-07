#pragma once
#include "initialweights.h"

class Customized : public InitialWeights
{
public:
    Customized(System *system);
    std::string getLabel() { return m_label; }
    void setupInitialWeights();

    Eigen::MatrixXd getParameters();

private:
    double m_factor = 1;
    std::string m_label = "customized";
};
