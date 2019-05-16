#include "rbmjastrow2.h"
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include <iostream>

RBMJastrow2::RBMJastrow2(System* system) :
        WaveFunction(system) {
    m_numberOfHiddenNodes               = m_system->getNumberOfHiddenNodes();
    m_numberOfFreeDimensions            = m_system->getNumberOfFreeDimensions();
    m_numberOfParameters                = m_numberOfParticles * m_numberOfParticles * m_numberOfHiddenNodes + m_numberOfHiddenNodes;
    double sigma                        = 1; //m_system->getWidth();
    m_sigmaSqrd                         = sigma*sigma;
    m_sigmaQuad                         = m_sigmaSqrd*m_sigmaSqrd;
}

void RBMJastrow2::setConstants(const int elementNumber) {
    m_maxNumberOfParametersPerElement   = m_system->getMaxNumberOfParametersPerElement();
    m_elementNumber                     = elementNumber;
}

void RBMJastrow2::updateVectors() {
    m_g = m_b1 + m_W1.transpose() * m_positions;
    m_h = m_b2 + m_W2.transpose() * m_positions + Eigen::VectorXd::Ones(m_numberOfHiddenNodes);
    Eigen::VectorXd m_e = m_g.cwiseProduct(m_h.cwiseInverse()).array().exp();
    m_p = (m_e + Eigen::VectorXd::Ones(m_numberOfHiddenNodes)).cwiseInverse();
    m_n = m_e.cwiseProduct(m_p);
}

void RBMJastrow2::updateRatio() {
    double Prod = 1;
    for(int j=0; j<m_numberOfHiddenNodes; j++) {
        Prod *= m_pOld(j)/m_p(j);
    }
    m_probabilityRatio  = Prod * Prod;
}

void RBMJastrow2::updateArrays(const Eigen::VectorXd positions, const Eigen::VectorXd radialVector, const Eigen::MatrixXd distanceMatrix, const int changedCoord) {
    m_positions = positions;
    updateVectors();
    updateRatio();
}

void RBMJastrow2::setArrays() {
    m_positionsOld          = m_positions;
    m_gOld                  = m_g;
    m_hOld                  = m_h;
    m_nOld                  = m_n;
    m_pOld                  = m_p;
    m_probabilityRatioOld   = m_probabilityRatio;
}

void RBMJastrow2::resetArrays() {
    m_positions             = m_positionsOld;
    m_h                     = m_hOld;
    m_g                     = m_gOld;
    m_n                     = m_nOld;
    m_p                     = m_pOld;
    m_probabilityRatio      = m_probabilityRatioOld;
}

void RBMJastrow2::initializeArrays(const Eigen::VectorXd positions, const Eigen::VectorXd radialVector, const Eigen::MatrixXd distanceMatrix) {
    m_positions         = positions;
    m_probabilityRatio  = 1;
    m_n  = Eigen::VectorXd::Zero(m_numberOfHiddenNodes);
    m_p  = Eigen::VectorXd::Zero(m_numberOfHiddenNodes);
    updateVectors();
}

void RBMJastrow2::updateParameters(Eigen::MatrixXd parameters) {
    m_b1     = parameters.row(m_elementNumber).head(m_numberOfHiddenNodes);
    m_b2     = parameters.row(m_elementNumber).segment(m_numberOfHiddenNodes,m_numberOfHiddenNodes);
    Eigen::VectorXd w1Flatten = parameters.row(m_elementNumber).segment(2*m_numberOfHiddenNodes, m_numberOfFreeDimensions*m_numberOfHiddenNodes);
    Eigen::Map<Eigen::MatrixXd> W1(w1Flatten.data(), m_numberOfFreeDimensions, m_numberOfHiddenNodes);
    m_W1     = W1;
    Eigen::VectorXd w2Flatten = parameters.row(m_elementNumber).segment(m_numberOfHiddenNodes*(2+m_numberOfFreeDimensions), m_numberOfFreeDimensions*m_numberOfHiddenNodes);
    Eigen::Map<Eigen::MatrixXd> W2(w2Flatten.data(), m_numberOfFreeDimensions, m_numberOfHiddenNodes);
    m_W2     = W2;
}

double RBMJastrow2::evaluateRatio() {
    return m_probabilityRatio;
}

double RBMJastrow2::computeGradient(const int k) {
    double Sum = 0;
    for(int j=0; j<m_numberOfHiddenNodes; j++) {
        Sum += m_n(j) * (m_W1(k,j) * m_h(j) - m_W2(k,j)*m_g(j)) / (m_h(j) * m_h(j));
    }
    return Sum / m_sigmaSqrd;
}

double RBMJastrow2::computeLaplacian() {
    double Sum = 0;
    for(int k=0; k<m_numberOfFreeDimensions; k++) {
        for(int j=0; j<m_numberOfHiddenNodes; j++) {
            Sum += m_n(j) * (m_p(j) * (m_W1(k,j) * m_h(j) - m_W2(k,j)*m_g(j)) * (m_W1(k,j) * m_h(j) - m_W2(k,j)*m_g(j)) / (m_h(j) * m_h(j) * m_h(j) * m_h(j))) - 2*m_W2(k,j) * (m_W1(k,j) * m_h(j) - m_W2(k,j)*m_g(j)) / (m_h(j) * m_h(j) * m_h(j));
        }
    }
    return Sum / m_sigmaQuad;
}

Eigen::VectorXd RBMJastrow2::computeParameterGradient() {
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);
    for(int l=0; l<m_numberOfHiddenNodes; l++) {
        gradients(l) = 1/m_h(l);
        gradients(l+m_numberOfHiddenNodes) = -1/(m_h(l) * m_h(l));
        for(int m=0; m<m_numberOfFreeDimensions; m++) {
            int n = l * m_numberOfFreeDimensions + m + 2*m_numberOfHiddenNodes;
            int o = l * m_numberOfFreeDimensions + m + m_numberOfHiddenNodes*(2 + m_numberOfFreeDimensions);
            gradients(n) = m_positions(m) / (m_h(l) * m_sigmaSqrd);
            gradients(o) = -m_positions(m) / (m_h(l) * m_h(l) * m_sigmaSqrd);
        }
    }
    return gradients;
}
