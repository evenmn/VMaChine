#include "randomnormal.h"
#include "../RNG/rng.h"
#include "../WaveFunctions/wavefunction.h"
#include "../system.h"

RandomNormal::RandomNormal(System *system)
    : InitialState(system)
{}

void RandomNormal::setupInitialState()
{
    m_numberOfDimensions = m_system->getNumberOfDimensions();
    m_numberOfParticles = m_system->getNumberOfParticles();
    m_degreesOfFreedom = m_system->getNumberOfFreeDimensions();
    m_omega = m_system->getFrequency();

    m_positions = m_system->getRandomNumberGenerator()->randomNormalMatrix(m_degreesOfFreedom, 1, 0, 1 / m_omega);
    InitialState::calculateDistanceMatrix();
    InitialState::calculateRadialVector();

    for (auto &i : m_system->getWaveFunctionElements()) {
        i->initializeArrays(m_positions, m_radialVector, m_distanceMatrix);
        i->setArrays();
    }
}
