#pragma once
#include <Eigen/Dense>
#include <iostream>
#include <cassert>

class InitialWeights
{
public:
    InitialWeights(class System *system);
    virtual std::string getLabel() = 0;
    virtual void setupInitialWeights() = 0;
    virtual Eigen::MatrixXd getParameters() = 0;

    virtual ~InitialWeights() = 0;

protected:
    int m_numberOfElements = 0;
    int m_maxParameters = 0;
    Eigen::MatrixXd m_parameters;

    class System *m_system = nullptr;
};
