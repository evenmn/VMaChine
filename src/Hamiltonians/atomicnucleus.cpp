#include "atomicnucleus.h"
#include "../WaveFunctions/wavefunction.h"
#include "../system.h"

AtomicNucleus::AtomicNucleus(System *system)
    : Hamiltonian(system)
{}

void AtomicNucleus::initialize()
{
    m_Z = m_system->getAtomicNumber();
    m_numberOfParticles = m_system->getNumberOfParticles();
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
