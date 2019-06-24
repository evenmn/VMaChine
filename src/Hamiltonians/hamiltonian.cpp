#include "hamiltonian.h"
#include "../system.h"

Hamiltonian::Hamiltonian(System* system) {
    m_system = system;
}

double Hamiltonian::getInteractionEnergy() {
    m_distanceMatrix        = m_system->getDistanceMatrix();
    double interactionEnergy = 0;
    if(m_interaction) {
        for(int i=0; i<m_numberOfParticles; i++) {
            for(int j=i+1; j<m_numberOfParticles; j++) {
                if(m_distanceMatrix(i,j) < m_dsl) {
                    interactionEnergy += 1/m_distanceMatrix(i,j);
                }
                else {
                    interactionEnergy += exp(-m_distanceMatrix(i,j)/m_screeningStrength)/m_distanceMatrix(i,j);
                }
            }
        }
    }
    return interactionEnergy;
}

Hamiltonian::~Hamiltonian() {}
