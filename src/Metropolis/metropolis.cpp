#include "metropolis.h"
#include "../system.h"

Metropolis::Metropolis(System *system)
{
    m_system = system;
    m_numberOfDimensions = m_system->getNumberOfDimensions();
    m_numberOfParticles = m_system->getNumberOfParticles();
    m_degreesOfFreedom = m_system->getNumberOfFreeDimensions();
    m_stepLength = m_system->getStepLength();
    m_positions = m_system->getInitialState()->getParticles();
    m_radialVector = m_system->getInitialState()->getRadialVector();
    m_distanceMatrix = m_system->getInitialState()->getDistanceMatrix();
    system->setGlobalArraysToCalculate();
    m_calculateDistanceMatrix = m_system->getCalculateDistanceMatrix();
    m_calculateRadialVector = m_system->getCalculateRadialVector();
    m_RNG = m_system->getRandomNumberGenerator();
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
    for (int i = 0; i < m_numberOfParticles; i++) {
        m_distanceMatrix(particle, i) = calculateDistanceMatrixElement(particle, i);
        m_distanceMatrix(i, particle) = m_distanceMatrix(particle, i);
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
