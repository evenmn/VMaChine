#include "doublewell.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../WaveFunctions/wavefunction.h"

HarmonicOscillator::HarmonicOscillator(System* system) :
        Hamiltonian(system) {
    m_omega                 = m_system->getFrequency();
    m_omega_sqrd            = m_omega * m_omega;
    assert(m_omega > 0);
    m_numberOfParticles     = m_system->getNumberOfParticles();
    m_numberOfDimensions    = m_system->getNumberOfDimensions();
    m_interaction           = m_system->getInteraction();
}

double HarmonicOscillator::computeLocalEnergy() {
    m_positions             = m_system->getPositions();
    m_distanceMatrix        = m_system->getDistanceMatrix();

    double interactionEnergy = 0;
    if(m_interaction) {
        for(int i=0; i<m_numberOfParticles; i++) {
            for(int j=i+1; j<m_numberOfParticles; j++) {
                interactionEnergy += 1/m_distanceMatrix(i,j);
            }
        }
    }

    double externalEnergy = 0.5 * m_omega_sqrd * (m_positions.cwiseAbs2()).sum();
    double kineticEnergy  = m_system->getKineticEnergy();

    return kineticEnergy + externalEnergy + interactionEnergy;
}
