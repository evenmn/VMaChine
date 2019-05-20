#include "partlyrestricted.h"
#include <cassert>
#include <iostream>
#include "../system.h"

PartlyRestricted::PartlyRestricted(System* system) :
        WaveFunction(system) {
    m_numberOfFreeDimensions            = m_system->getNumberOfFreeDimensions();
    m_numberOfParameters                = m_numberOfFreeDimensions * m_numberOfFreeDimensions;
}

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

void PartlyRestricted::setConstants(const unsigned int elementNumber) {
    m_maxNumberOfParametersPerElement   = m_system->getMaxParameters();
    m_elementNumber                     = elementNumber;
}

void PartlyRestricted::calculateProbabilityRatio(int changedCoord) {
    //double ratio = (m_beta.row(particle).transpose() * (m_h.row(particle) - m_hOld.row(particle))).sum();

    double ratio = 0;
    for(int i=0; i<changedCoord+1; i++) {
        ratio -= abs(m_positions(i) * m_c(i,changedCoord) * m_positions(changedCoord));
    }
    for(int j=changedCoord; j<m_numberOfFreeDimensions; j++) {
        ratio -= abs(m_positions(changedCoord) * m_c(changedCoord,j) * m_positions(j));
    }
    m_probabilityRatio = exp(ratio);
}

void PartlyRestricted::updateArrays(const Eigen::VectorXd positions, const Eigen::VectorXd radialVector, const Eigen::MatrixXd distanceMatrix, const int changedCoord) {
    m_positions = positions;
    //m_xCx = positions.transpose() * m_c * positions;
    calculateProbabilityRatio(changedCoord);
}

void PartlyRestricted::setArrays() {
    m_positionsOld          = m_positions;
    //m_xCxOld                = m_xCx;
    m_probabilityRatioOld   = m_probabilityRatio;
}

void PartlyRestricted::resetArrays() {
    m_positions             = m_positionsOld;
    //m_xCx                   = m_xCxOld;
    m_probabilityRatio      = m_probabilityRatioOld;
}

void PartlyRestricted::initializeArrays(const Eigen::VectorXd positions, const Eigen::VectorXd radialVector, const Eigen::MatrixXd distanceMatrix) {
    m_positions = positions;
    //m_xCx = positions.transpose() * m_c * positions;

    m_probabilityRatio = 1;
}

void PartlyRestricted::updateParameters(Eigen::MatrixXd parameters) {
    Eigen::Map<Eigen::MatrixXd> c(parameters.row(m_elementNumber).data(), m_numberOfFreeDimensions, m_numberOfFreeDimensions);
    m_c = c;
    //m_c = Eigen::MatrixXd::Zero(m_numberOfFreeDimensions, m_numberOfFreeDimensions);
}

double PartlyRestricted::evaluateRatio() {
    return m_probabilityRatio;
}

double PartlyRestricted::computeGradient(const int k) {
    double Sum = 0;
    for(int i=0; i<k+1; i++) {
        Sum -= m_positions(i) * m_c(i,k) * sgn(m_positions(i) * m_c(i,k) * m_positions(k));
    }
    for(int j=k; j<m_numberOfFreeDimensions; j++) {
        Sum -= m_c(k,j) * m_positions(j) * sgn(m_positions(k) * m_c(k,j) * m_positions(j));
    }
    return Sum;
}

double PartlyRestricted::computeLaplacian() {
    double Sum = 0;
    for(int i=0; i<m_numberOfFreeDimensions; i++) {
        Sum -= m_c(i,i) * sgn(m_positions(i) * m_c(i,i) * m_positions(i));
    }
    return 2 * Sum;
}

Eigen::VectorXd PartlyRestricted::computeParameterGradient() {
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);

    for(int m=0; m<m_numberOfFreeDimensions; m++) {
        for(int l=m; l<m_numberOfFreeDimensions; l++) {
            gradients(l * m_numberOfFreeDimensions + m) = -m_positions(m) * m_positions(l) * sgn(m_positions(m) * m_c(m,l) * m_positions(l));
        }
    }
    return gradients;
}
