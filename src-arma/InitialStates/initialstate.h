#pragma once
#include <armadillo>

class InitialState
{
public:
    InitialState(class System *system);
    virtual void setupInitialState() = 0;
    arma::vec getParticles() { return m_positions; }
    arma::mat getDistanceMatrix() { return m_distanceMatrix; }
    arma::vec getRadialVector() { return m_radialVector; }

    double calculateDistanceMatrixElement(arma::uword i, arma::uword j);
    void calculateDistanceMatrix();
    double calculateRadialVectorElement(arma::uword particle);
    void calculateRadialVector();

    virtual ~InitialState() = 0;

protected:
    class System *m_system = nullptr;
    arma::uword m_degreesOfFreedom = 0;
    arma::uword m_numberOfDimensions = 0;
    arma::uword m_numberOfParticles = 0;

    arma::vec m_positions;
    arma::mat m_distanceMatrix;
    arma::vec m_radialVector;
};
