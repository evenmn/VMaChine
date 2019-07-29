#include "atomicnucleus.h"
#include "../WaveFunctions/wavefunction.h"
#include "../system.h"

AtomicNucleus::AtomicNucleus(System *system)
    : Hamiltonian(system)
{
    m_Z = m_system->getAtomicNumber();
    m_numberOfParticles = m_system->getNumberOfParticles();
    m_numberOfDimensions = m_system->getNumberOfDimensions();
    m_interaction = m_system->getInteraction();
    m_screening = m_system->getScreening();
    m_screeningStrength = m_system->getScreeningStrength();
    m_dsl = m_system->getDSL();
}

double AtomicNucleus::getExternalEnergy()
{
    m_positions = m_system->getPositions();
    m_radialVector = m_system->getRadialVector();
    double nucleusEnergy = 0;
    for (int i = 0; i < m_numberOfParticles; i++) {
        nucleusEnergy += 1 / m_radialVector(i);
    }
    int l = 0;
    return -m_Z * nucleusEnergy + l * (l + 1) * m_positions.cwiseAbs2().cwiseInverse().sum() / 2;
}

double AtomicNucleus::computeLocalEnergy()
{
    double kineticEnergy = m_system->getKineticEnergy();
    double externalEnergy = getExternalEnergy();
    double interactionEnergy = Hamiltonian::getInteractionEnergy();
    return kineticEnergy + externalEnergy + interactionEnergy;
}
