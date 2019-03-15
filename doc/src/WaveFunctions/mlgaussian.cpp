#include "mlgaussian.h"
#include <cassert>
#include <iostream>
#include "../system.h"

MLGaussian::MLGaussian(System* system) :
        WaveFunction(system) {
    m_numberOfFreeDimensions            = m_system->getNumberOfFreeDimensions();
    m_maxNumberOfParametersPerElement   = m_system->getMaxNumberOfParametersPerElement();
    m_omega                             = m_system->getFrequency();
    double sigma                        = m_system->getWidth();
    m_sigmaSqrd = sigma*sigma;
}

void MLGaussian::updateParameters(const Eigen::MatrixXd parameters, const int elementNumber) {
    m_elementNumber = elementNumber;
    m_a = (parameters.row(m_elementNumber)).head(m_numberOfFreeDimensions);
}

void MLGaussian::initializeArrays(const Eigen::VectorXd positions) {
    m_positions = positions;
    m_Xa        = positions - m_a;
    m_ratio     = 1;
}

void MLGaussian::updateArrays(const Eigen::VectorXd positions, const int pRand) {
    m_oldPositions  = m_positions;
    m_positions     = positions;
    m_oldXa         = m_Xa;
    m_Xa            = positions - m_a;
    m_oldRatio      = m_ratio;
    m_ratio         = exp((m_oldXa(pRand)*m_oldXa(pRand) - m_Xa(pRand)*m_Xa(pRand)) / (2 * m_sigmaSqrd));
}

void MLGaussian::resetArrays(int pRand) {
    m_positions = m_oldPositions;
    m_Xa        = m_oldXa;
    m_ratio     = m_oldRatio;
}

double MLGaussian::evaluateRatio() {
    return m_ratio;
}

double MLGaussian::computeFirstDerivative(const int k) {
    return - m_omega * m_Xa(k)/m_sigmaSqrd;
}

double MLGaussian::computeSecondDerivative() {;
    return - m_omega * m_numberOfFreeDimensions/m_sigmaSqrd;
}

Eigen::VectorXd MLGaussian::computeFirstEnergyDerivative(const int k) {
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);
    gradients(k) = - 0.5 * m_omega / m_sigmaSqrd;
    return gradients;
}

Eigen::VectorXd MLGaussian::computeSecondEnergyDerivative() {
    return Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);
}
