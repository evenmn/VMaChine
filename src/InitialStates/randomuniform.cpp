#include "randomuniform.h"
#include "../RNG/rng.h"
#include "../WaveFunctions/wavefunction.h"
#include "../system.h"

RandomUniform::RandomUniform(System *system)
    : InitialState(system)
{}

void RandomUniform::setupInitialState()
{
    m_numberOfParticles = m_system->getNumberOfParticles();
    m_numberOfDimensions = m_system->getNumberOfDimensions();
    m_degreesOfFreedom = m_system->getNumberOfFreeDimensions();

    m_positions = m_system->getRandomNumberGenerator()->randomUniformMatrix(m_degreesOfFreedom, 1);
    InitialState::calculateDistanceMatrix();
    InitialState::calculateRadialVector();

    for (auto &i : m_system->getWaveFunctionElements()) {
        i->initializeArrays(m_positions, m_radialVector, m_distanceMatrix);
        i->setArrays();
    }
}
