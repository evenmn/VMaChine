#pragma once
#include <armadillo>
#include <iostream>
#include <string>

class WaveFunction
{
public:
    WaveFunction(class System *system);
    virtual arma::uword getNumberOfParameters() = 0;
    virtual arma::uword getGlobalArrayNeed() = 0;
    virtual std::string getLabel() = 0;

    virtual void updateParameters(const arma::mat parameters) = 0;
    virtual void initializeArrays(const arma::vec positions,
                                  const arma::vec radialVector,
                                  const arma::mat distanceMatrix)
        = 0;
    virtual void updateArrays(const arma::vec positions,
                              const arma::vec radialVector,
                              const arma::mat distanceMatrix,
                              const arma::uword changedCoord)
        = 0;
    virtual void setConstants(const arma::uword elementNumber) = 0;
    virtual void setArrays() = 0;
    virtual void resetArrays() = 0;
    virtual double evaluateRatio() = 0;
    virtual double computeGradient(const arma::uword k) = 0;
    virtual double computeLaplacian() = 0;
    virtual arma::vec computeParameterGradient() = 0;

    virtual ~WaveFunction() = 0;

protected:
    arma::uword m_numberOfParticles = 0;
    arma::uword m_numberOfDimensions = 0;
    arma::uword m_degreesOfFreedom = 0;

    class System *m_system = nullptr;
};
