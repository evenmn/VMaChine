#include "atomicnucleus.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;

AtomicNucleus::AtomicNucleus(System* system) :
        Hamiltonian(system) {
    m_Z                     = m_system->getAtomicNumber();
    m_numberOfParticles     = m_system->getNumberOfParticles();
    m_numberOfDimensions    = m_system->getNumberOfDimensions();
    m_interaction           = m_system->getInteraction();
}

double AtomicNucleus::computeLocalEnergy() {
    m_positions             = m_system->getParticles();

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

    double nucleusEnergy = 0;
    for(int i=0; i<m_numberOfParticles; i++) {
        double sqrdElementWise = 0;
        for(int d=0; d<m_numberOfDimensions; d++) {
            sqrdElementWise += m_positions(i*m_numberOfDimensions + d) * m_positions(i*m_numberOfDimensions + d);
        }
        nucleusEnergy += 1/sqrt(sqrdElementWise);
    }

    int l = 0;
    double externalEnergy = - m_Z * nucleusEnergy + l*(l+1) * m_positions.cwiseAbs2().cwiseInverse().sum()/2;
    double kineticEnergy  = m_system->getKineticEnergy();

    return kineticEnergy + externalEnergy + interactionEnergy;
}
