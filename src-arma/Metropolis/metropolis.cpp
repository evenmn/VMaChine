#include "metropolis.h"
#include "../system.h"

Metropolis::Metropolis(System *system)
{
    m_system = system;
}

double Metropolis::calculateDistanceMatrixElement(const arma::uword i, const arma::uword j)
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

void Metropolis::calculateDistanceMatrixCross(const arma::uword particle)
{
    for (arma::uword i = 0; i < m_numberOfParticles; i++) {
        m_distanceMatrix(particle, i) = calculateDistanceMatrixElement(particle, i);
        m_distanceMatrix(i, particle) = m_distanceMatrix(particle, i);
    }
}

void Metropolis::calculateRadialVectorElement(const arma::uword particle)
{
    double sqrdElementWise = 0;
    arma::uword initPartCoord = particle * m_numberOfDimensions;
    for (arma::uword i = initPartCoord; i < initPartCoord + m_numberOfDimensions; i++) {
        sqrdElementWise += m_positions(i) * m_positions(i);
    }
    m_radialVector(particle) = sqrt(sqrdElementWise);
}

Metropolis::~Metropolis() {}
