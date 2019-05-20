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
    m_numberOfElements    = m_system->getNumberOfElements();
    m_maxParameters = m_system->getMaxParameters();
    m_waveFunctionVector              = m_system->getWaveFunctionElements();
    m_eta                             = m_system->getLearningRate();
    m_parameters                      = m_system->getInitialWeights()->getWeights();
    m_gradients                       = Eigen::MatrixXd::Zero(m_numberOfElements, m_maxParameters);
    m_oldParameters                   = Eigen::MatrixXd::Zero(m_numberOfElements, m_maxParameters);
}

Eigen::MatrixXd BarzilaiBorwein::updateParameters() {
    m_oldGradients = m_gradients;
    m_parameters    = m_system->getWeights();
    m_gradients     = Optimization::getEnergyGradient();
    Eigen::MatrixXd learningRate = (m_parameters - m_oldParameters).cwiseProduct((m_gradients - m_oldGradients).cwiseInverse());
    for(int i=0; i<m_numberOfElements; i++) {
        for(int j=0; j<m_maxParameters; j++) {
            if(std::isnan(learningRate(i,j)) || std::isinf(learningRate(i,j))) {
                learningRate(i,j) = m_eta;
            }
        }
    }
    m_oldParameters = m_parameters;
    return m_eta * learningRate.cwiseProduct(m_gradients);
}
