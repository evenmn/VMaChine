#include "doublewell.h"
#include "../WaveFunctions/wavefunction.h"
#include "../system.h"
#include <iostream>

DoubleWell::DoubleWell(System *system, double displacement)
    : Hamiltonian(system)
{
    m_b = displacement;
}

void DoubleWell::initialize()
{
    double omega = m_system->getFrequency();
    m_omega_sqrd = omega * omega;
    m_numberOfParticles = m_system->getNumberOfParticles();
    m_numberOfDimensions = m_system->getNumberOfDimensions();
    m_offset = m_numberOfParticles * m_b * m_b / 4;
}

double DoubleWell::getExternalEnergy()
{
    m_positions = m_system->getPositions();
    double sumX = 0;
    for (int i = 0; i < m_numberOfParticles; i++) {
        sumX += fabs(m_positions(i * m_numberOfDimensions));
    }

    return 0.5 * m_omega_sqrd * (m_positions.cwiseAbs2().sum() - m_b * fabs(sumX) + m_offset);
}
