#include "initialstate.h"

InitialState::InitialState(System *system)
{
    m_system = system;
}

double InitialState::calculateDistanceMatrixElement(arma::uword i, arma::uword j)
{
    double dist = 0;
    arma::uword parti = m_numberOfDimensions * i;
    arma::uword partj = m_numberOfDimensions * j;
    for (arma::uword d = 0; d < m_numberOfDimensions; d++) {
        double diff = m_positions(parti + d) - m_positions(partj + d);
        dist += diff * diff;
    }
    return sqrt(dist);
}

void InitialState::calculateDistanceMatrix()
{
    m_distanceMatrix.zeros(m_numberOfParticles, m_numberOfParticles);
    for (arma::uword i = 0; i < m_numberOfParticles; i++) {
        for (arma::uword j = i + 1; j < m_numberOfParticles; j++) {
            m_distanceMatrix(i, j) = calculateDistanceMatrixElement(i, j);
            m_distanceMatrix(j, i) = m_distanceMatrix(i, j);
        }
    }
}

double InitialState::calculateRadialVectorElement(arma::uword particle)
{
    double sqrtElementWise = 0;
    arma::uword part = particle * m_numberOfDimensions;
    for (arma::uword d = 0; d < m_numberOfDimensions; d++) {
        sqrtElementWise += m_positions(part + d) * m_positions(part + d);
    }
    return sqrt(sqrtElementWise);
}

void InitialState::calculateRadialVector()
{
    m_radialVector.zeros(m_numberOfParticles);
    for (arma::uword i = 0; i < m_numberOfParticles; i++) {
        m_radialVector(i) = calculateRadialVectorElement(i);
    }
}

InitialState::~InitialState() {}
