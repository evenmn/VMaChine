#include "bruteforce.h"
#include "../InitialStates/initialstate.h"
#include "../RNG/rng.h"
#include "../WaveFunctions/wavefunction.h"
#include "../system.h"
#include <cassert>
#include <iostream>

BruteForce::BruteForce(System *system)
    : Metropolis(system)
{
    m_numberOfDimensions = m_system->getNumberOfDimensions();
    m_numberOfParticles = m_system->getNumberOfParticles();
    m_numberOfFreeDimensions = m_system->getNumberOfFreeDimensions();
    m_stepLength = m_system->getStepLength();
    m_positions = m_system->getInitialState()->getParticles();
    m_radialVector = m_system->getInitialState()->getRadialVector();
    m_distanceMatrix = m_system->getInitialState()->getDistanceMatrix();
    system->setGlobalArraysToCalculate();
    m_calculateDistanceMatrix = m_system->getCalculateDistanceMatrix();
    m_calculateRadialVector = m_system->getCalculateRadialVector();
}

bool BruteForce::acceptMove()
{
    int pRand = m_system->getRandomNumberGenerator()->nextInt(m_numberOfFreeDimensions);

    m_positionsOld = m_positions;
    m_radialVectorOld = m_radialVector;
    m_distanceMatrixOld = m_distanceMatrix;

    m_positions(pRand) = m_positionsOld(pRand)
                         + 10 * (m_system->getRandomNumberGenerator()->nextDouble() - 0.5)
                               * m_stepLength;
    if (m_calculateDistanceMatrix) {
        Metropolis::calculateDistanceMatrixCross(int(pRand / m_numberOfDimensions));
    }
    if (m_calculateRadialVector) {
        Metropolis::calculateRadialVectorElement(int(pRand / m_numberOfDimensions));
    }
    m_system->updateAllArrays(m_positions, m_radialVector, m_distanceMatrix, pRand);

    double ratio = m_system->evaluateWaveFunctionRatio();
    double r = m_system->getRandomNumberGenerator()->nextDouble();
    if (ratio < r) {
        m_positions(pRand) = m_positionsOld(pRand);
        m_distanceMatrix = m_distanceMatrixOld;
        m_radialVector = m_radialVectorOld;
        m_system->resetAllArrays();
        return false;
    }
    return true;
}
