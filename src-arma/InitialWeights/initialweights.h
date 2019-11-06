#pragma once
#include <armadillo>
#include <iostream>

class InitialWeights
{
public:
    InitialWeights(class System *system);
    virtual void setupInitialWeights() = 0;
    virtual arma::mat getParameters() = 0;

    virtual ~InitialWeights() = 0;

protected:
    arma::uword m_numberOfElements = 0;
    arma::uword m_maxParameters = 0;
    arma::mat m_parameters;

    class System *m_system = nullptr;
};
