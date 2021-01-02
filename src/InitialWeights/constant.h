#pragma once
#include "initialweights.h"

class Constant : public InitialWeights
{
public:
    Constant(System *system, const double factor = 1);
    std::string getLabel() { return m_label; }
    void setupInitialWeights();

    Eigen::MatrixXd getParameters();

private:
    double m_factor = 1;
    std::string m_label = "constant";
};
