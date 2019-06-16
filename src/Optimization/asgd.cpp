#include "asgd.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../sampler.h"
#include "../WaveFunctions/wavefunction.h"

ASGD::ASGD(System* system, const double gamma) :
        Optimization(system) {
    m_numberOfFreeDimensions    = m_system->getNumberOfFreeDimensions();
    m_numberOfElements          = m_system->getNumberOfElements();
    m_maxParameters             = m_system->getMaxParameters();
    m_waveFunctionVector        = m_system->getWaveFunctionElements();
    m_eta                       = m_system->getLearningRate();
    m_v                         = Eigen::MatrixXd::Ones(m_numberOfElements, m_maxParameters);
    m_gamma                     = gamma;
}

Eigen::MatrixXd ASGD::updateParameters() {
    m_step += 1;
    //m_t = m_t + m_A;
    m_v = m_gamma * m_v + m_eta * Optimization::getEnergyGradient() / sqrt(m_step);
    return m_v;
}
