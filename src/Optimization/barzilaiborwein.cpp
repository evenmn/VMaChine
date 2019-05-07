#include "barzilaiborwein.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../sampler.h"
#include "../WaveFunctions/wavefunction.h"
#include "../InitialWeights/initialweights.h"

BarzilaiBorwein::BarzilaiBorwein(System* system) :
        Optimization(system) {
    m_numberOfFreeDimensions          = m_system->getNumberOfFreeDimensions();
    m_numberOfWaveFunctionElements    = m_system->getNumberOfWaveFunctionElements();
    m_maxNumberOfParametersPerElement = m_system->getMaxNumberOfParametersPerElement();
    m_waveFunctionVector              = m_system->getWaveFunctionElements();
    m_eta                             = m_system->getLearningRate();
    m_parameters                      = m_system->getInitialWeights()->getWeights();
    m_gradients                       = Eigen::MatrixXd::Zero(m_numberOfWaveFunctionElements, m_maxNumberOfParametersPerElement);
    m_oldParameters                   = Eigen::MatrixXd::Zero(m_numberOfWaveFunctionElements, m_maxNumberOfParametersPerElement);
}

Eigen::MatrixXd BarzilaiBorwein::getEnergyGradient() {
    double          averageEnergy     = m_system->getSampler()->getAverageEnergy();
    Eigen::MatrixXd averageGradients  = m_system->getSampler()->getAverageGradients();
    Eigen::MatrixXd averageGradientsE = m_system->getSampler()->getAverageGradientsE();
    return 2 * (averageGradientsE - averageEnergy * averageGradients);
}

Eigen::MatrixXd BarzilaiBorwein::updateParameters() {
    m_oldGradients  = m_gradients;
    m_parameters    = m_system->getWeights();
    m_gradients     = getEnergyGradient();
    Eigen::MatrixXd learningRate = (m_parameters - m_oldParameters).cwiseProduct((m_gradients - m_oldGradients).cwiseInverse());
    for(unsigned int i=0; i<m_numberOfWaveFunctionElements; i++) {
        for(unsigned int j=0; j<m_maxNumberOfParametersPerElement; j++) {
            if(std::isnan(learningRate(i,j)) || std::isinf(learningRate(i,j))) {
                learningRate(i,j) = m_eta;
            }
        }
    }
    m_oldParameters = m_parameters;
    return learningRate.cwiseProduct(m_gradients) * Optimization::getEnergyGradient();
}
