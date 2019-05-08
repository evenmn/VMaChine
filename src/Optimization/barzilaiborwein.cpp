#include "barzilaiborwein.h"
#include "../system.h"
#include "../InitialWeights/initialweights.h"

BarzilaiBorwein::BarzilaiBorwein(System* system) :
        Optimization(system) {
    m_numberOfWaveFunctionElements    = m_system->getNumberOfWaveFunctionElements();
    m_maxNumberOfParametersPerElement = m_system->getMaxNumberOfParametersPerElement();
    m_eta                             = m_system->getLearningRate();
    m_parameters                      = m_system->getInitialWeights()->getWeights();
    m_gradients                       = Eigen::MatrixXd::Zero(m_numberOfWaveFunctionElements, m_maxNumberOfParametersPerElement);
    m_oldParameters                   = Eigen::MatrixXd::Zero(m_numberOfWaveFunctionElements, m_maxNumberOfParametersPerElement);
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
