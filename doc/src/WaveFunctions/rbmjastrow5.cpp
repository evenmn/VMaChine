#include "rbmjastrow5.h"
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include <iostream>

RBMJastrow5::RBMJastrow5(System* system) :
        WaveFunction(system) {
    m_numberOfHiddenNodes               = m_system->getNumberOfHiddenNodes();
    m_numberOfFreeDimensions            = m_system->getNumberOfFreeDimensions();
    m_maxNumberOfParametersPerElement   = m_system->getMaxNumberOfParametersPerElement();
    double sigma                        = 1; //m_system->getWidth();
    m_sigmaSqrd                         = sigma*sigma;
    m_sigmaQuad                         = m_sigmaSqrd*m_sigmaSqrd;
}

void RBMJastrow5::updateGradient() {
    for(int k=0; k<m_numberOfFreeDimensions; k++) {
        for(int j=0; j<m_numberOfHiddenNodes; j++) {
            m_gradientPart(k,j) = m_W1(k,j) / m_sigmaSqrd + 2 * m_positions(k) * m_W2(k,j) / m_sigmaQuad;
        }
    }
}

void RBMJastrow5::updateVectors() {
    m_v = 2*m_b + m_W1.transpose() * m_positions/m_sigmaSqrd + m_W2.transpose() * m_positionsSqrd/m_sigmaQuad;
    Eigen::VectorXd m_e = m_v.array().exp();
    m_p = (m_e + Eigen::VectorXd::Ones(m_numberOfHiddenNodes)).cwiseInverse();
    m_n = m_e.cwiseProduct(m_p);

    m_pDotN   = m_p.cwiseProduct(m_n);
}

void RBMJastrow5::updateRatio() {
    double Prod = 1;
    for(int j=0; j<m_numberOfHiddenNodes; j++) {
        Prod *= m_pOld(j)/m_p(j);
    }
    m_probabilityRatio  = Prod * Prod;
}

void RBMJastrow5::updateArrays(const Eigen::VectorXd positions, const int changedCoord) {
    setArrays();
    m_positions     = positions;
    m_positionsSqrd = positions.cwiseAbs2();
    updateVectors();
    updateRatio();
    updateGradient();
}

void RBMJastrow5::setArrays() {
    m_positionsOld          = m_positions;
    m_positionsSqrdOld      = m_positionsSqrd;
    m_gradientPartOld       = m_gradientPart;
    m_vOld                  = m_v;
    m_nOld                  = m_n;
    m_pOld                  = m_p;
    m_pDotNOld              = m_pDotN;
    m_pMinusNOld            = m_pMinusN;
    m_probabilityRatioOld   = m_probabilityRatio;
}

void RBMJastrow5::resetArrays() {
    m_positions             = m_positionsOld;
    m_positionsSqrd         = m_positionsSqrdOld;
    m_gradientPart          = m_gradientPartOld;
    m_v                     = m_vOld;
    m_n                     = m_nOld;
    m_p                     = m_pOld;
    m_pDotN                 = m_pDotNOld;
    m_pMinusN               = m_pMinusNOld;
    m_probabilityRatio      = m_probabilityRatioOld;
}

void RBMJastrow5::initializeArrays(const Eigen::VectorXd positions) {
    m_positions         = positions;
    m_positionsSqrd     = positions.cwiseAbs2();
    m_probabilityRatio  = 1;

    m_n             = Eigen::VectorXd::Zero(m_numberOfHiddenNodes);
    m_p             = Eigen::VectorXd::Zero(m_numberOfHiddenNodes);
    m_gradientPart  = Eigen::MatrixXd::Zero(m_numberOfFreeDimensions, m_numberOfHiddenNodes);

    updateVectors();
    updateGradient();
    setArrays();
}

void RBMJastrow5::updateParameters(Eigen::MatrixXd parameters, const int elementNumber) {
    m_elementNumber = elementNumber;
    Eigen::VectorXd w1Flatten = parameters.row(m_elementNumber).segment(m_numberOfHiddenNodes, m_numberOfFreeDimensions*m_numberOfHiddenNodes);
    Eigen::Map<Eigen::MatrixXd> W1(w1Flatten.data(), m_numberOfFreeDimensions, m_numberOfHiddenNodes);
    Eigen::VectorXd w2Flatten = parameters.row(m_elementNumber).segment(m_numberOfHiddenNodes*(1+m_numberOfFreeDimensions), m_numberOfFreeDimensions*m_numberOfHiddenNodes);
    Eigen::Map<Eigen::MatrixXd> W2(w2Flatten.data(), m_numberOfFreeDimensions, m_numberOfHiddenNodes);
    m_W1     = W1;
    m_W2     = W2;
    m_b     = parameters.row(m_elementNumber).head(m_numberOfHiddenNodes);
}

double RBMJastrow5::evaluateRatio() {
    return m_probabilityRatio;
}

double RBMJastrow5::computeGradient(const int k) {
    double sum = 0;
    for(int j=0; j<m_numberOfHiddenNodes; j++) {
        sum += m_n(j) * m_gradientPart(k,j);
    }
    return sum;
}

double RBMJastrow5::computeLaplacian() {
    double sum = 0;
    for(int k=0; k<m_numberOfFreeDimensions; k++) {
        for(int j=0; j<m_numberOfHiddenNodes; j++) {
            sum += m_n(j) * (2*m_W2(k,j) / m_sigmaQuad + m_p(j) * m_gradientPart(k,j) * m_gradientPart(k,j));
        }
    }
    return sum;
}

Eigen::VectorXd RBMJastrow5::computeParameterGradient() {
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);
    for(int l=0; l<m_numberOfHiddenNodes; l++) {
        gradients(l) = 2 * m_n(l);
        for(int m=0; m<m_numberOfFreeDimensions; m++) {
            int n = l * m_numberOfFreeDimensions + m + m_numberOfHiddenNodes;
            int o = l * m_numberOfFreeDimensions + m + m_numberOfHiddenNodes * (1 + m_numberOfFreeDimensions);
            gradients(n) = m_positions(m) * m_n(l) / m_sigmaSqrd;
            gradients(o) = m_positionsSqrd(m) * m_n(l) / m_sigmaSqrd;
        }
    }
    return gradients;
}
