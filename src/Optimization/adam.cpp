#include "adam.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../sampler.h"
#include "../WaveFunctions/wavefunction.h"

ADAM::ADAM(System* system) :
        Optimization(system) {
    m_numberOfFreeDimensions          = m_system->getNumberOfFreeDimensions();
    m_numberOfWaveFunctionElements    = m_system->getNumberOfWaveFunctionElements();
    m_maxNumberOfParametersPerElement = m_system->getMaxNumberOfParametersPerElement();
    m_waveFunctionVector              = m_system->getWaveFunctionElements();
    m_eta                             = m_system->getLearningRate();
    m_m                               = Eigen::MatrixXd::Ones(m_numberOfWaveFunctionElements, m_maxNumberOfParametersPerElement);
    m_v                               = Eigen::MatrixXd::Ones(m_numberOfWaveFunctionElements, m_maxNumberOfParametersPerElement);
    m_theta                           = Eigen::MatrixXd::Ones(m_numberOfWaveFunctionElements, m_maxNumberOfParametersPerElement);
}

Eigen::MatrixXd ADAM::updateParameters() {
    m_step += 1;
    m_g = Optimization::getEnergyGradient();
    m_m = m_beta1 * m_m + (1 - m_beta1) * m_g;
    m_v = m_beta2 * m_v + (1 - m_beta2) * m_g.cwiseAbs2();
    m_mHat = m_m/(1 - pow(m_beta1, m_step));
    m_vHat = m_v/(1 - pow(m_beta2, m_step));
    for(unsigned int i=0; i<m_numberOfWaveFunctionElements; i++) {
        for(unsigned int j=0; j<m_maxNumberOfParametersPerElement; j++) {
            m_theta(i,j) = m_eta * m_mHat(i,j)/(sqrt(m_vHat(i,j) + m_epsilon));
        }
    }
    return m_theta;
}
