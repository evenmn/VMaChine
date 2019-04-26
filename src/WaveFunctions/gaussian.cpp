#include "gaussian.h"
#include <cassert>
#include <iostream>
#include "../system.h"

Gaussian::Gaussian(System* system) :
        WaveFunction(system) {
    m_numberOfFreeDimensions            = m_system->getNumberOfFreeDimensions();
    m_maxNumberOfParametersPerElement   = m_system->getMaxNumberOfParametersPerElement();
    m_omega                             = m_system->getFrequency();
}

void Gaussian::initializeArrays(const Eigen::VectorXd positions) {
    m_positions             = positions;
    m_probabilityRatio      = 1;
    setArrays();
}

void Gaussian::updateProbabilityRatio(int changedCoord) {
    m_probabilityRatio      = exp(m_omega * m_alpha * (m_positionsOld(changedCoord)*m_positionsOld(changedCoord) \
                                                     - m_positions   (changedCoord)*m_positions   (changedCoord)));
}


void Gaussian::updateArrays(const Eigen::VectorXd positions, const int changedCoord) {
    setArrays();
    m_positions             = positions;
    updateProbabilityRatio(changedCoord);
}

void Gaussian::setArrays() {
    m_positionsOld          = m_positions;
    m_probabilityRatioOld   = m_probabilityRatio;
}

void Gaussian::resetArrays() {
    m_positions             = m_positionsOld;
    m_probabilityRatio      = m_probabilityRatioOld;
}

void Gaussian::updateParameters(const Eigen::MatrixXd parameters, const int elementNumber) {
    m_elementNumber         = elementNumber;
    m_alpha                 = parameters(m_elementNumber,0);
}

double Gaussian::evaluateRatio() {
    return m_probabilityRatio;
}

double Gaussian::computeGradient(const int k) {
    return - m_omega * m_alpha * m_positions(k);
}

double Gaussian::computeLaplacian() {;
    return - m_omega * m_alpha * m_numberOfFreeDimensions;
}

Eigen::VectorXd Gaussian::computeParameterGradient() {
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);
    gradients(0) = - 0.5 * m_omega * m_positions.cwiseAbs2().sum();
    return gradients;
}
