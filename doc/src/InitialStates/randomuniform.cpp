#include "randomuniform.h"
#include <iostream>
#include <cassert>
#include "RNG/rng.h"
#include "WaveFunctions/wavefunction.h"
#include "../system.h"

RandomUniform::RandomUniform(System*    system)  :
        InitialState(system) {
    m_numberOfFreeDimensions = m_system->getNumberOfFreeDimensions();
    setupInitialState();
}

void RandomUniform::setupInitialState() {
    m_positions = Eigen::VectorXd::Zero(m_numberOfFreeDimensions);
    for (int i=0; i < m_numberOfFreeDimensions; i++) {
        m_positions(i) = m_system->getRandomNumberGenerator()->nextDouble();
    }

    for(auto& i : m_system->getWaveFunctionElements()) {
        i->initializeArrays(m_positions);
    }
}
