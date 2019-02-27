#include "bruteforce.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../WaveFunctions/wavefunction.h"
#include "../RNG/rng.h"
#include "InitialStates/initialstate.h"

BruteForce::BruteForce(System* system) :
        Metropolis(system) {
    m_numberOfFreeDimensions = m_system->getNumberOfFreeDimensions();
    m_stepLength             = m_system->getStepLength();
    m_positions              = m_system->getInitialState()->getParticles();
}

bool BruteForce::acceptMove() {
    double psiOld       = m_system->evaluateWaveFunctionSqrd();
    int pRand = m_system->getRandomNumberGenerator()->nextInt(m_numberOfFreeDimensions);

    m_positionsOld      = m_positions;
    m_positions(pRand) = m_positionsOld(pRand) + 10 * (m_system->getRandomNumberGenerator()->nextDouble() - 0.5) * m_stepLength;

    m_system->updateAllArrays(m_positions, pRand);
    double psiNew = m_system->evaluateWaveFunctionSqrd();

    double w = psiNew/psiOld;
    double r = m_system->getRandomNumberGenerator()->nextDouble();
    if(w < r) {
        m_positions(pRand) = m_positionsOld(pRand);
        m_system->resetAllArrays();
        return false;
    }
    return true;
}
