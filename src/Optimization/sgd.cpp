#include "sgd.h"
#include <iostream>
#include "../system.h"
#include "../WaveFunctions/wavefunction.h"

SGD::SGD(System* system, const double gamma, const double monotonicExp) :
        Optimization(system) {
    m_numberOfWaveFunctionElements    = m_system->getNumberOfWaveFunctionElements();
    m_maxNumberOfParametersPerElement = m_system->getMaxNumberOfParametersPerElement();
    m_eta                             = m_system->getLearningRate();
    m_v                               = Eigen::MatrixXd::Ones(m_numberOfWaveFunctionElements, m_maxNumberOfParametersPerElement);
    m_gamma                           = gamma;
    m_monotonicExp                    = monotonicExp;
}

Eigen::MatrixXd SGD::updateParameters() {
    m_iter += 1;
    double monotonic = 1/pow(m_iter, m_monotonicExp);
    m_v = m_gamma * m_v + m_eta * Optimization::getEnergyGradient() * monotonic;
    return m_v;
}
