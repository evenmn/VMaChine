#pragma once
#include "initialweights.h"

class Javier: public InitialWeights
{
public:
    Javier(System *system);
    std::string getLabel() { return m_label; }
    void setupInitialWeights();

    Eigen::MatrixXd getParameters();

private:
    std::string m_label = "javier";
};
