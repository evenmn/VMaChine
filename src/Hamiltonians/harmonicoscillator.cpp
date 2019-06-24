#include "harmonicoscillator.h"
#include <cassert>
#include "../system.h"
#include "../WaveFunctions/wavefunction.h"

HarmonicOscillator::HarmonicOscillator(System* system) :
        Hamiltonian(system) {
    m_omega                 = m_system->getFrequency();
    m_omegaSqrd             = m_omega * m_omega;
    m_numberOfParticles     = m_system->getNumberOfParticles();
    m_numberOfDimensions    = m_system->getNumberOfDimensions();
    m_interaction           = m_system->getInteraction();
    m_screeningStrength     = m_system->getScreeningStrength();
    m_dsl                   = m_system->getDSL();
}

double HarmonicOscillator::getExternalEnergy() {
    m_positions = m_system->getPositions();
    return 0.5 * m_omegaSqrd * m_positions.cwiseAbs2().sum();
}

double HarmonicOscillator::computeLocalEnergy() {
    double kineticEnergy     = m_system->getKineticEnergy();
    double externalEnergy    = getExternalEnergy();
    double interactionEnergy = Hamiltonian::getInteractionEnergy();
    return kineticEnergy + externalEnergy + interactionEnergy;
}

