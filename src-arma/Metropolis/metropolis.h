#pragma once
#include <armadillo>

class Metropolis
{
public:
    Metropolis(class System *system);
    virtual arma::vec updatePositions() { return m_positions; }
    virtual arma::vec updateRadialVector() { return m_radialVector; }
    virtual arma::mat updateDistanceMatrix() { return m_distanceMatrix; }

    virtual void initialize() = 0;
    virtual bool acceptMove() = 0;
    virtual ~Metropolis() = 0;

    double calculateDistanceMatrixElement(const arma::uword i, const arma::uword j);
    void calculateDistanceMatrixCross(const arma::uword particle);
    void calculateRadialVectorElement(const arma::uword particle);

protected:
    class System *m_system = nullptr;
    arma::vec m_positions;
    arma::vec m_positionsOld;
    arma::vec m_radialVector;
    arma::vec m_radialVectorOld;
    arma::mat m_distanceMatrix;
    arma::mat m_distanceMatrixOld;
    arma::uword m_degreesOfFreedom = 0;
    arma::uword m_numberOfDimensions = 0;
    arma::uword m_numberOfParticles = 0;
    double m_stepLength = 0;
    double m_diff = 0.5;

    bool m_calculateDistanceMatrix = false;
    bool m_calculateRadialVector = false;

    class RandomNumberGenerator *m_RNG = nullptr;
};
