#include "partlyrestricted.h"
#include <cassert>
#include <iostream>
#include "../system.h"

PartlyRestricted::PartlyRestricted(System* system) :
        WaveFunction(system) {
    m_numberOfFreeDimensions            = m_system->getNumberOfFreeDimensions();
    m_maxNumberOfParametersPerElement   = m_system->getMaxNumberOfParametersPerElement();
    double sigma                        = m_system->getWidth();
    m_sigmaSqrd2 = sigma*sigma*sigma*sigma;
}

void PartlyRestricted::updateArrays(const Eigen::VectorXd positions, const int pRand) {
    m_oldPositions = m_positions;
    m_positions = positions;

    m_oldXCx = m_xCx;
    m_xCx = positions.transpose() * m_c * positions;
}

void PartlyRestricted::resetArrays() {
    m_positions = m_oldPositions;
    m_xCx       = m_oldXCx;
}

void PartlyRestricted::initializeArrays(const Eigen::VectorXd positions) {
    m_positions = positions;
    m_xCx = positions.transpose() * m_c * positions;
}

void PartlyRestricted::updateParameters(Eigen::MatrixXd parameters, const int elementNumber) {
    m_elementNumber = elementNumber;
    Eigen::Map<Eigen::MatrixXd> c(parameters.row(m_elementNumber).data(), m_numberOfFreeDimensions, m_numberOfFreeDimensions);
    m_c = c;
}

double PartlyRestricted::evaluate() {
    //return exp(-0.5 * m_xCx  / m_sigmaSqrd2);

    double Sum = 0;
    for(int i=0; i<m_numberOfFreeDimensions; i++) {
        for(int j=i; j<m_numberOfFreeDimensions; j++) {
            Sum += m_positions(i) * m_positions(j) * m_c(i,j);
        }
    }
    return exp(-0.5 * Sum);
}

double PartlyRestricted::evaluateSqrd() {
    //return exp(- m_xCx  / m_sigmaSqrd2);

    double Sum = 0;
    for(int i=0; i<m_numberOfFreeDimensions; i++) {
        for(int j=i; j<m_numberOfFreeDimensions; j++) {
            Sum += m_positions(i) * m_positions(j) * m_c(i,j);
        }
    }
    return exp(-Sum);
}

double PartlyRestricted::computeFirstDerivative(const int k) {
    //return - double(m_c.row(k) * positions) / m_sigmaSqrd2;

    double Sum = 0;
    for(int j=k+1; j<m_numberOfFreeDimensions; j++) {
        Sum += m_positions(j) * m_c(k,j);
    }
    return -0.5 * Sum - m_positions(k) * m_c(k,k);
}

double PartlyRestricted::computeSecondDerivative() {
    return -m_c.diagonal().sum();
}

Eigen::VectorXd PartlyRestricted::computeFirstEnergyDerivative(const int k) {
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);

    for(int j=0; j<m_numberOfFreeDimensions; j++) {
        if(j==k) {
            gradients(j * m_numberOfFreeDimensions + k) = 0.5 * m_positions(k);
        }
        else {
            gradients(j * m_numberOfFreeDimensions + k) = 0.25 * m_positions(j);
        }
    }
    return gradients;
}

Eigen::VectorXd PartlyRestricted::computeSecondEnergyDerivative() {
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);

    for(int j=0; j<m_numberOfFreeDimensions; j++) {
        gradients(j * m_numberOfFreeDimensions + j) = 0.5;
    }
    return gradients;
}
