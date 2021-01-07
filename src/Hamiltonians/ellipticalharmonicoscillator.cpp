#include "ellipticalharmonicoscillator.h"
#include "../WaveFunctions/wavefunction.h"
#include "../system.h"
#include <cassert>
#include <cmath>

EllipticalHarmonicOscillator::EllipticalHarmonicOscillator(System *system, const double beta)
    : Hamiltonian(system)
{
    m_beta = beta;
}

void EllipticalHarmonicOscillator::initialize()
{
    double omega = m_system->getFrequency();
    m_omegaSqrd = omega * omega;
}

double EllipticalHarmonicOscillator::getExternalEnergy()
{
    m_positions = m_system->getPositions();
    Eigen::MatrixXd positionsSqrd = m_positions.cwiseAbs2();
    double sum = 0;
    for (int i=0; i<m_positions.size(); i++) {
        if (i+1 % 3 == 0) {
            sum += m_beta * positionsSqrd(i);
        } else {
            sum += positionsSqrd(i);
        }
    }
    return 0.5 * m_omegaSqrd * sum;
}
