#pragma once
#include <string>

#include "initialweights.h"

class FromFile : public InitialWeights
{
public:
    FromFile(System *system, std::string = "weights.dat");
    std::string getLabel() { return m_label; }
    void setupInitialWeights();

    Eigen::MatrixXd getParameters();

private:
    double m_factor = 1;
    std::string m_filename;
    std::string m_label = "from-file";
};
