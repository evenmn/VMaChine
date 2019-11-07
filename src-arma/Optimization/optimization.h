#pragma once
#include <armadillo>
#include <vector>

class Optimization
{
public:
    Optimization(class System *system);
    virtual arma::uword getNumberOfBatches() = 0;
    virtual std::string getLabel() = 0;

    virtual void initialize() = 0;
    virtual arma::mat updateParameters() = 0;
    virtual ~Optimization() = 0;

    arma::mat getEnergyGradient();

protected:
    class System *m_system = nullptr;
    arma::vec m_positions;
    arma::uword m_numberOfElements = 0;
    arma::uword m_maxParameters = 0;
    double m_eta = 1;
    std::vector<class WaveFunction *> m_waveFunctionVector;
};
