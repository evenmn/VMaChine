#include "doublewell.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../WaveFunctions/wavefunction.h"

DoubleWell::DoubleWell(System* system) :
        Hamiltonian(system) {
    m_omega                 = m_system->getFrequency();
    m_omega_sqrd            = m_omega * m_omega;
    assert(m_omega > 0);
    m_numberOfParticles     = m_system->getNumberOfParticles();
    m_numberOfDimensions    = m_system->getNumberOfDimensions();
    m_interaction           = m_system->getInteraction();
}

double DoubleWell::computeLocalEnergy() {
    m_positions             = m_system->getPositions();
    m_distanceMatrix        = m_system->getDistanceMatrix();

    double interactionEnergy = 0;
    if(m_interaction) {
        for(unsigned int i=0; i<m_numberOfParticles; i++) {
            for(unsigned int j=i+1; j<m_numberOfParticles; j++) {
                interactionEnergy += 1/m_distanceMatrix(i,j);
            }
        }
    }

    double sumX = 0;
    for(unsigned int i=0; i<m_numberOfParticles; i++) {
        sumX += m_positions(int(i*m_numberOfDimensions));
    }

    double externalEnergy = 0;
    if(m_numberOfDimensions == 2) {
        externalEnergy = 0.5 * m_omega_sqrd * ((m_positions.cwiseAbs2()).sum() + 0.25 * m_A * m_A - m_A * sumX);
    }
    else {
        std::cout << "DoubleWell potential is implemented in two dimensions only" << std::endl;
        MPI_Finalize();
        exit(0);
    }
    double kineticEnergy  = m_system->getKineticEnergy();

    return kineticEnergy + externalEnergy + interactionEnergy;
}
