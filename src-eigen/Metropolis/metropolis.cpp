#include "metropolis.h"
#include "../system.h"

Metropolis::Metropolis(System *system)
{
    m_system = system;
}

double Metropolis::calculateDistanceMatrixElement(const int i, const int j)
{
    double dist = 0;
    int parti = m_numberOfDimensions * i;
    int partj = m_numberOfDimensions * j;
    for (int d = 0; d < m_numberOfDimensions; d++) {
        double diff = m_positions(parti + d) - m_positions(partj + d);
        dist += diff * diff;
    }
    return sqrt(dist);
}

void Metropolis::calculateDistanceMatrixCross(const int particle)
{
    for (int i = 0; i < particle; i++) {
        m_distanceMatrix(i, particle) = calculateDistanceMatrixElement(particle, i);
    }
    for (int i = particle+1; i < m_numberOfParticles; i++) {
        m_distanceMatrix(particle, i) = calculateDistanceMatrixElement(particle, i);
    }
}

void Metropolis::calculateRadialVectorElement(const int particle)
{
    double sqrdElementWise = 0;
    int initPartCoord = particle * m_numberOfDimensions;
    for (int i = initPartCoord; i < initPartCoord + m_numberOfDimensions; i++) {
        sqrdElementWise += m_positions(i) * m_positions(i);
    }
    m_radialVector(particle) = sqrt(sqrdElementWise);
}

Metropolis::~Metropolis() {}
