#include "harmonicoscillator.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;

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

    double interactionEnergy = 0;
    if(m_interaction) {
        for(int i=0; i<m_numberOfParticles; i++) {
            for(int j=0; j<i; j++) {
                double sqrdElementWise = 0;
                for(int d=0; d<m_numberOfDimensions; d++) {
                    double numb = m_positions(i*m_numberOfDimensions + d) - m_positions(j*m_numberOfDimensions + d);
                    sqrdElementWise += numb * numb;
                }
                interactionEnergy += 1/sqrt(sqrdElementWise);
            }
        }
    }

    double externalEnergy = 0.5 * m_omega_sqrd * (m_positions.cwiseAbs2()).sum();
    double kineticEnergy  = m_system->getKineticEnergy();

    return kineticEnergy + externalEnergy + interactionEnergy;
}

