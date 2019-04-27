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
    m_positions             = m_system->getPositions();
    m_distanceMatrix        = m_system->getDistanceMatrix();
    m_radialVector          = m_system->getRadialVector();

    double interactionEnergy = 0;
    if(m_interaction) {
        for(int i=0; i<m_numberOfParticles; i++) {
            for(int j=i+1; j<m_numberOfParticles; j++) {
                interactionEnergy += 1/m_distanceMatrix(i,j);
            }
        }
    }

    double nucleusEnergy = 0;
    for(int i=0; i<m_numberOfParticles; i++) {
        nucleusEnergy += 1/m_radialVector(i);
    }

    int l = 0;
    double externalEnergy = - m_Z * nucleusEnergy + l*(l+1) * m_positions.cwiseAbs2().cwiseInverse().sum()/2;
    double kineticEnergy  = m_system->getKineticEnergy();

    return kineticEnergy + externalEnergy + interactionEnergy;
}
