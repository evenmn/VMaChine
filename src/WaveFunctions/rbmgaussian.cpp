#include "rbmgaussian.h"
#include <cassert>
#include <iostream>
#include "../system.h"

RBMGaussian::RBMGaussian(System* system) :
        WaveFunction(system) {
    m_numberOfFreeDimensions            = m_system->getNumberOfFreeDimensions();
    m_maxNumberOfParametersPerElement   = m_system->getMaxNumberOfParametersPerElement();
    m_omega                             = m_system->getFrequency();
    double sigma                        = m_system->getWidth();
    m_sigmaSqrd = sigma*sigma;
}

void RBMGaussian::updateParameters(const Eigen::MatrixXd parameters, const int elementNumber) {
    m_elementNumber         = elementNumber;
    m_a                     = (parameters.row(m_elementNumber)).head(m_numberOfFreeDimensions);
}

void RBMGaussian::initializeArrays(const Eigen::VectorXd positions) {
    m_positions             = positions;
    m_Xa                    = positions - m_a;
    m_probabilityRatio      = 1;

    setArrays();
}

void RBMGaussian::updateArrays(const Eigen::VectorXd positions, const int changedCoord) {
    setArrays();

    m_positions             = positions;
    m_Xa                    = positions - m_a;
    m_probabilityRatio      = exp((m_XaOld(changedCoord)*m_XaOld(changedCoord) - m_Xa(changedCoord)*m_Xa(changedCoord)) / (2 * m_sigmaSqrd));
}

void RBMGaussian::setArrays() {
    m_positionsOld          = m_positions;
    m_XaOld                 = m_Xa;
    m_probabilityRatioOld   = m_probabilityRatio;
}

void RBMGaussian::resetArrays() {
    m_positions             = m_positionsOld;
    m_Xa                    = m_XaOld;
    m_probabilityRatio      = m_probabilityRatioOld;
}

double RBMGaussian::evaluateRatio() {
    return m_probabilityRatio;
}

double RBMGaussian::computeGradient(const int k) {
    return - m_Xa(k)/m_sigmaSqrd;
}

double RBMGaussian::computeLaplacian() {;
    return - m_numberOfFreeDimensions/m_sigmaSqrd;
}

Eigen::VectorXd RBMGaussian::computeParameterGradient() {
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);
    gradients.head(m_numberOfFreeDimensions) = m_Xa/m_sigmaSqrd;
    return gradients;
}