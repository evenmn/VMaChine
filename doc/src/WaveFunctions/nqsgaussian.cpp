#include "nqsgaussian.h"
#include <cassert>
#include <iostream>
#include "../system.h"

NQSGaussian::NQSGaussian(System* system) :
        WaveFunction(system) {
    m_numberOfFreeDimensions            = m_system->getNumberOfFreeDimensions();
    m_maxNumberOfParametersPerElement   = m_system->getMaxNumberOfParametersPerElement();
    m_omega                             = m_system->getFrequency();
    double sigma                        = m_system->getWidth();
    m_sigmaSqrd = sigma*sigma;
}

void NQSGaussian::updateParameters(const Eigen::MatrixXd parameters, const int elementNumber) {
    m_elementNumber         = elementNumber;
    m_a                     = (parameters.row(m_elementNumber)).head(m_numberOfFreeDimensions);
}

void NQSGaussian::initializeArrays(const Eigen::VectorXd positions) {
    m_positions             = positions;
    m_Xa                    = positions - m_a;
    m_probabilityRatio      = 1;
}

void NQSGaussian::updateArrays(const Eigen::VectorXd positions, const int pRand) {
    m_positionsOld          = m_positions;
    m_positions             = positions;
    m_XaOld                 = m_Xa;
    m_Xa                    = positions - m_a;
    m_probabilityRatioOld   = m_probabilityRatio;
    m_probabilityRatio      = exp((m_XaOld(pRand)*m_XaOld(pRand) - m_Xa(pRand)*m_Xa(pRand)) / (2 * m_sigmaSqrd));
}

void NQSGaussian::resetArrays(int pRand) {
    m_positions             = m_positionsOld;
    m_Xa                    = m_XaOld;
    m_probabilityRatio      = m_probabilityRatioOld;
}

double NQSGaussian::evaluateRatio() {
    return m_probabilityRatio;
}

double NQSGaussian::computeFirstDerivative(const int k) {
    return - m_Xa(k)/m_sigmaSqrd;
}

double NQSGaussian::computeSecondDerivative() {;
    return - m_numberOfFreeDimensions/m_sigmaSqrd;
}

Eigen::VectorXd NQSGaussian::computeFirstEnergyDerivative(const int k) {
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);
    gradients(k) = - 0.5 / m_sigmaSqrd;
    return gradients;
}

Eigen::VectorXd NQSGaussian::computeSecondEnergyDerivative() {
    return Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);
}