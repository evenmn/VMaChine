#include "randomuniform.h"
#include "../RNG/rng.h"
#include "../WaveFunctions/wavefunction.h"
#include "../system.h"
#include <cassert>
#include <iostream>

RandomUniform::RandomUniform(System *system)
    : InitialState(system)
{
    m_numberOfParticles = m_system->getNumberOfParticles();
    m_numberOfDimensions = m_system->getNumberOfDimensions();
    m_degreesOfFreedom = m_system->getNumberOfFreeDimensions();
    setupInitialState();
}

double RandomUniform::calculateDistanceMatrixElement(int i, int j)
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

void RandomUniform::calculateDistanceMatrix()
{
    m_distanceMatrix = Eigen::MatrixXd::Zero(m_numberOfParticles, m_numberOfParticles);
    for (int i = 0; i < m_numberOfParticles; i++) {
        for (int j = i + 1; j < m_numberOfParticles; j++) {
            m_distanceMatrix(i, j) = calculateDistanceMatrixElement(i, j);
            m_distanceMatrix(j, i) = m_distanceMatrix(i, j);
        }
    }
}

double RandomUniform::calculateRadialVectorElement(int particle)
{
    double sqrtElementWise = 0;
    int part = particle * m_numberOfDimensions;
    for (int d = 0; d < m_numberOfDimensions; d++) {
        sqrtElementWise += m_positions(part + d) * m_positions(part + d);
    }
    return sqrt(sqrtElementWise);
}

void RandomUniform::calculateRadialVector()
{
    m_radialVector = Eigen::VectorXd::Zero(m_numberOfParticles);
    for (int i = 0; i < m_numberOfParticles; i++) {
        m_radialVector(i) = calculateRadialVectorElement(i);
    }
}

void RandomUniform::setupInitialState()
{
    m_positions = Eigen::VectorXd::Zero(m_degreesOfFreedom);
    for (int i = 0; i < m_degreesOfFreedom; i++) {
        m_positions(i) = m_system->getRandomNumberGenerator()->nextDouble();
    }
    InitialState::calculateDistanceMatrix();
    InitialState::calculateRadialVector();

    for (auto &i : m_system->getWaveFunctionElements()) {
        i->initializeArrays(m_positions, m_radialVector, m_distanceMatrix);
        i->setArrays();
    }
}
