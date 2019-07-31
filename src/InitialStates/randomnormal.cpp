#include "randomnormal.h"
#include "../RNG/rng.h"
#include "../WaveFunctions/wavefunction.h"
#include "../system.h"
#include <cassert>
#include <iostream>

RandomNormal::RandomNormal(System *system)
    : InitialState(system)
{
    m_numberOfParticles = m_system->getNumberOfParticles();
    m_numberOfDimensions = m_system->getNumberOfDimensions();
    m_degreesOfFreedom = m_system->getNumberOfFreeDimensions();
    setupInitialState();
}

double RandomNormal::calculateDistanceMatrixElement(int i, int j)
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

void RandomNormal::calculateDistanceMatrix()
{
    m_distanceMatrix = Eigen::MatrixXd::Zero(m_numberOfParticles, m_numberOfParticles);
    for (int i = 0; i < m_numberOfParticles; i++) {
        for (int j = i + 1; j < m_numberOfParticles; j++) {
            m_distanceMatrix(i, j) = calculateDistanceMatrixElement(i, j);
            m_distanceMatrix(j, i) = m_distanceMatrix(i, j);
        }
    }
}

double RandomNormal::calculateRadialVectorElement(int particle)
{
    double sqrtElementWise = 0;
    int part = particle * m_numberOfDimensions;
    for (int d = 0; d < m_numberOfDimensions; d++) {
        sqrtElementWise += m_positions(part + d) * m_positions(part + d);
    }
    return sqrt(sqrtElementWise);
}

void RandomNormal::calculateRadialVector()
{
    m_radialVector = Eigen::VectorXd::Zero(m_numberOfParticles);
    for (int i = 0; i < m_numberOfParticles; i++) {
        m_radialVector(i) = calculateRadialVectorElement(i);
    }
}

void RandomNormal::setupInitialState()
{
    m_positions = Eigen::VectorXd::Zero(m_degreesOfFreedom);
    for (int i = 0; i < m_degreesOfFreedom; i++) {
        m_positions(i) = m_system->getRandomNumberGenerator()->nextGaussian(0, 1);
    }
    InitialState::calculateDistanceMatrix();
    InitialState::calculateRadialVector();
    for (auto &i : m_system->getWaveFunctionElements()) {
        i->initializeArrays(m_positions, m_radialVector, m_distanceMatrix);
        i->setArrays();
    }
}
