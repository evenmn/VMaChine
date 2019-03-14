#include "asgd.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../sampler.h"
#include "../WaveFunctions/wavefunction.h"

ASGD::ASGD(System* system, const double gamma) :
        Optimization(system) {
    m_numberOfFreeDimensions          = m_system->getNumberOfFreeDimensions();
    m_numberOfWaveFunctionElements    = m_system->getNumberOfWaveFunctionElements();
    m_maxNumberOfParametersPerElement = m_system->getMaxNumberOfParametersPerElement();
    m_waveFunctionVector              = m_system->getWaveFunctionElements();
    m_eta                             = m_system->getLearningRate();
    m_v                               = Eigen::MatrixXd::Ones(m_numberOfWaveFunctionElements, m_maxNumberOfParametersPerElement);
    m_gamma                           = gamma;
}

Eigen::VectorXd ASGD::getInstantGradients(WaveFunction* waveFunction) {
    Eigen::VectorXd TotalGradients = waveFunction->computeSecondEnergyDerivative();
    for(int k = 0; k < m_numberOfFreeDimensions; k++) {
        double Sum = 0;
        for(auto& i : m_waveFunctionVector) {
            Sum += i->computeFirstDerivative(k);
        }
        TotalGradients += 2 * Sum * waveFunction->computeFirstEnergyDerivative(k);
    }
    return TotalGradients;
}

Eigen::MatrixXd ASGD::getAllInstantGradients() {
    Eigen::MatrixXd gradients = Eigen::MatrixXd::Zero(m_numberOfWaveFunctionElements, m_maxNumberOfParametersPerElement);
    for(int i = 0; i < m_numberOfWaveFunctionElements; i++) {
        gradients.row(i) += getInstantGradients(m_waveFunctionVector[unsigned(i)]);
    }
    return gradients;
}

Eigen::MatrixXd ASGD::getEnergyGradient() {
    double          averageEnergy     = m_system->getSampler()->getAverageEnergy();
    Eigen::MatrixXd averageGradients  = m_system->getSampler()->getAverageGradients();
    Eigen::MatrixXd averageGradientsE = m_system->getSampler()->getAverageGradientsE();
    return 2 * (averageGradientsE - averageEnergy * averageGradients);
}

Eigen::MatrixXd ASGD::updateParameters() {
    m_step += 1;
    //m_t = m_t + m_A;
    m_v = m_gamma * m_v + m_eta * getEnergyGradient() / sqrt(m_step);
    return m_v;
}
