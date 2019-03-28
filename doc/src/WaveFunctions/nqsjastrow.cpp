#include "nqsjastrow.h"
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include <iostream>

NQSJastrow::NQSJastrow(System* system) :
        WaveFunction(system) {
    m_numberOfHiddenNodes               = m_system->getNumberOfHiddenNodes();
    m_numberOfFreeDimensions            = m_system->getNumberOfFreeDimensions();
    m_maxNumberOfParametersPerElement   = m_system->getMaxNumberOfParametersPerElement();
    double sigma                        = 1; //m_system->getWidth();
    m_sigmaSqrd                         = sigma*sigma;
    m_sigmaQuad                         = m_sigmaSqrd*m_sigmaSqrd;
}

void NQSJastrow::updateVectors() {
    m_v = m_b + m_W.transpose() * m_positions;
    Eigen::VectorXd m_e = m_v.array().exp();
    m_p = (m_e + Eigen::VectorXd::Ones(m_numberOfHiddenNodes)).cwiseInverse();
    m_n = m_e.cwiseProduct(m_p);

    m_pDotN   = m_p.cwiseProduct(m_n);
    m_pMinusN = m_p - m_n;
}

void NQSJastrow::updateRatio() {
    double Prod = 1;
    for(int j=0; j<m_numberOfHiddenNodes; j++) {
        Prod *= m_pOld(j)/m_p(j);
    }
    m_probabilityRatio  = Prod * Prod;
}

void NQSJastrow::updateArrays(const Eigen::VectorXd positions, const int pRand) {
    setArrays();
    m_positions = positions;
    updateVectors();
    updateRatio();
}

void NQSJastrow::setArrays() {
    m_positionsOld          = m_positions;
    m_vOld                  = m_v;
    m_nOld                  = m_n;
    m_pOld                  = m_p;
    m_pDotNOld              = m_pDotN;
    m_pMinusNOld            = m_pMinusN;
    m_probabilityRatioOld   = m_probabilityRatio;
}

void NQSJastrow::resetArrays() {
    m_positions             = m_positionsOld;
    m_v                     = m_vOld;
    m_n                     = m_nOld;
    m_p                     = m_pOld;
    m_pDotN                 = m_pDotNOld;
    m_pMinusN               = m_pMinusNOld;
    m_probabilityRatio      = m_probabilityRatioOld;
}

void NQSJastrow::initializeArrays(const Eigen::VectorXd positions) {
    m_positions         = positions;
    m_probabilityRatio  = 1;

    m_n  = Eigen::VectorXd::Zero(m_numberOfHiddenNodes);
    m_p  = Eigen::VectorXd::Zero(m_numberOfHiddenNodes);

    updateVectors();
    setArrays();
}

void NQSJastrow::updateParameters(Eigen::MatrixXd parameters, const int elementNumber) {
    m_elementNumber = elementNumber;
    Eigen::VectorXd wFlatten = parameters.row(m_elementNumber).segment(m_numberOfHiddenNodes, m_numberOfFreeDimensions*m_numberOfHiddenNodes);
    Eigen::Map<Eigen::MatrixXd> W(wFlatten.data(), m_numberOfFreeDimensions, m_numberOfHiddenNodes);
    m_W     = W;
    m_WSqrd = W.cwiseAbs2();
    m_b     = parameters.row(m_elementNumber).head(m_numberOfHiddenNodes);
}

double NQSJastrow::evaluateRatio() {
    return m_probabilityRatio;
}

double NQSJastrow::computeGradient(const int k) {
    //return double(m_W.row(k) * m_n)/sigma_sqrd;

    double Sum = 0;
    for(int j=0; j<m_numberOfHiddenNodes; j++) {
        Sum += m_W(k,j) * m_n(j);
    }
    return Sum / m_sigmaSqrd;
}

double NQSJastrow::computeLaplacian() {
    //return (m_W.cwiseAbs2()*m_p.cwiseProduct(m_n)).sum();

    double Sum = 0;
    for(int k=0; k<m_numberOfFreeDimensions; k++) {
        for(int j=0; j<m_numberOfHiddenNodes; j++) {
            Sum += m_WSqrd(k,j) * m_pDotN(j);
        }
    }
    return Sum / m_sigmaQuad;
}

Eigen::VectorXd NQSJastrow::computeParameterGradient() {
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);
    for(int l=0; l<m_numberOfHiddenNodes; l++) {
        gradients(l) = m_n(l);
        for(int m=0; m<m_numberOfFreeDimensions; m++) {
            int n = l * m_numberOfFreeDimensions + m + m_numberOfHiddenNodes;
            gradients(n) = m_positions(m) * m_n(l) / m_sigmaSqrd;
        }
    }
    return gradients;
}
