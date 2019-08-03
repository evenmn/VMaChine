#include "bruteforce.h"
#include "../InitialStates/initialstate.h"
#include "../RNG/rng.h"
#include "../WaveFunctions/wavefunction.h"
#include "../system.h"
#include <cassert>
#include <iostream>

BruteForce::BruteForce(System *system)
    : Metropolis(system)
{}

bool BruteForce::acceptMove()
{
    int i = m_RNG->nextInt(m_degreesOfFreedom);

    m_positionsOld = m_positions;
    m_radialVectorOld = m_radialVector;
    m_distanceMatrixOld = m_distanceMatrix;

    m_positions(i) += 10 * (m_RNG->nextDouble() - 0.5) * m_stepLength;
    if (m_calculateDistanceMatrix) {
        Metropolis::calculateDistanceMatrixCross(int(i / m_numberOfDimensions));
    }
    if (m_calculateRadialVector) {
        Metropolis::calculateRadialVectorElement(int(i / m_numberOfDimensions));
    }
    m_system->updateAllArrays(m_positions, m_radialVector, m_distanceMatrix, i);

    double p = m_system->evaluateProbabilityRatio();
    if (p < m_RNG->nextDouble()) {
        m_positions = m_positionsOld;
        m_distanceMatrix = m_distanceMatrixOld;
        m_radialVector = m_radialVectorOld;
        m_system->resetAllArrays();
        return false;
    }
    return true;
}
