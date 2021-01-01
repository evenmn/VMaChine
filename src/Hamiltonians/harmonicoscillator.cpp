#include "harmonicoscillator.h"
#include "../WaveFunctions/wavefunction.h"
#include "../system.h"
#include <cassert>

HarmonicOscillator::HarmonicOscillator(System *system)
    : Hamiltonian(system)
{}

void HarmonicOscillator::initialize()
{
    double omega = m_system->getFrequency();
    m_omegaSqrd = omega * omega;
}

double HarmonicOscillator::getExternalEnergy()
{
    m_positions = m_system->getPositions();
    return 0.5 * m_omegaSqrd * m_positions.cwiseAbs2().sum();
}
