#include "hamiltonian.h"
#include "../system.h"

Hamiltonian::Hamiltonian(System *system)
{
    m_system = system;
}

double Hamiltonian::getInteractionEnergy()
{
    m_distanceMatrix = m_system->getDistanceMatrix();
    double interactionEnergy = 0;
    if (m_interaction) {
        //interactionEnergy = m_distanceMatrix.cwiseInverse().sum() / 2;

        for (arma::uword i = 0; i < m_numberOfParticles; i++) {
            for (arma::uword j = i + 1; j < m_numberOfParticles; j++) {
                if (m_distanceMatrix(i, j) > m_dsl && m_screening) {
                    interactionEnergy += exp(-m_distanceMatrix(i, j) / m_screeningStrength)
                                         / m_distanceMatrix(i, j);
                } else {
                    interactionEnergy += 1 / m_distanceMatrix(i, j);
                }
            }
        }

    }
    return interactionEnergy;
}

Hamiltonian::~Hamiltonian() {}
