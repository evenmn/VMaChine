#include "asgd.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../sampler.h"
#include "../WaveFunctions/wavefunction.h"

ASGD::ASGD(System* system, const double gamma) :
        Optimization(system) {
    m_numberOfWaveFunctionElements    = m_system->getNumberOfWaveFunctionElements();
    m_maxNumberOfParametersPerElement = m_system->getMaxNumberOfParametersPerElement();
    m_eta                             = m_system->getLearningRate();
    m_v                               = Eigen::MatrixXd::Ones(m_numberOfWaveFunctionElements, m_maxNumberOfParametersPerElement);
    m_gamma                           = gamma;
}

Eigen::MatrixXd ASGD::updateParameters() {
    m_iter += 1;
    //m_t = m_t + m_A;
    m_v = m_gamma * m_v + m_eta * Optimization::getEnergyGradient() / sqrt(m_iter);
    return m_v;
}
