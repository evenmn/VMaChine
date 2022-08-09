#pragma once
#include "initialweights.h"

class Xavier: public InitialWeights
{
public:
    Xavier(System *system);
    std::string getLabel() { return m_label; }
    void setupInitialWeights();

    Eigen::MatrixXd getParameters();

private:
    std::string m_label = "Xavier";
};
