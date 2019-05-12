#include "rbmjastrow.h"
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include <iostream>

RBMJastrow::RBMJastrow(System* system) :
        WaveFunction(system) {
    m_numberOfHiddenNodes               = m_system->getNumberOfHiddenNodes();
    m_numberOfFreeDimensions            = m_system->getNumberOfFreeDimensions();
    m_numberOfParameters                = m_numberOfHiddenNodes * m_numberOfFreeDimensions + m_numberOfHiddenNodes;
    double sigma                        = 1; //m_system->getWidth();
    m_sigmaSqrd                         = sigma*sigma;
    m_sigmaQuad                         = m_sigmaSqrd*m_sigmaSqrd;
}

void RBMJastrow::updateParameters(Eigen::MatrixXd parameters, const int elementNumber) {
    m_elementNumber                     = elementNumber;
    m_maxNumberOfParametersPerElement   = m_system->getMaxNumberOfParametersPerElement();
    Eigen::VectorXd wFlatten = parameters.row(m_elementNumber).segment(m_numberOfHiddenNodes, m_numberOfFreeDimensions*m_numberOfHiddenNodes);
    Eigen::Map<Eigen::MatrixXd> W(wFlatten.data(), m_numberOfFreeDimensions, m_numberOfHiddenNodes);
    m_W     = W;
    m_WSqrd = W.cwiseAbs2();
    m_b     = parameters.row(m_elementNumber).head(m_numberOfHiddenNodes);
}

void RBMJastrow::initializeArrays(const Eigen::VectorXd positions, const Eigen::VectorXd radialVector, const Eigen::MatrixXd distanceMatrix) {
    m_positions         = positions;
    m_probabilityRatio  = 1;
    m_n  = Eigen::VectorXd::Zero(m_numberOfHiddenNodes);
    m_p  = Eigen::VectorXd::Zero(m_numberOfHiddenNodes);
    updateVectors();
}

void RBMJastrow::updateArrays(const Eigen::VectorXd positions, const Eigen::VectorXd radialVector, const Eigen::MatrixXd distanceMatrix, const int changedCoord) {
    m_positions = positions;
    updateVectors();
    updateRatio();
}

void RBMJastrow::setArrays() {
    m_positionsOld          = m_positions;
    m_vOld                  = m_v;
    m_nOld                  = m_n;
    m_pOld                  = m_p;
    m_pDotNOld              = m_pDotN;
    m_probabilityRatioOld   = m_probabilityRatio;
}

void RBMJastrow::resetArrays() {
    m_positions             = m_positionsOld;
    m_v                     = m_vOld;
    m_n                     = m_nOld;
    m_p                     = m_pOld;
    m_pDotN                 = m_pDotNOld;
    m_probabilityRatio      = m_probabilityRatioOld;
}

double RBMJastrow::evaluateRatio() {
    return m_probabilityRatio;
}

double RBMJastrow::computeGradient(const int k) {
    double Sum = 0;
    for(int j=0; j<m_numberOfHiddenNodes; j++) {
        Sum += m_W(k,j) * m_n(j);
    }
    return Sum / m_sigmaSqrd;
}

double RBMJastrow::computeLaplacian() {
    double Sum = 0;
    for(int k=0; k<m_numberOfFreeDimensions; k++) {
        for(int j=0; j<m_numberOfHiddenNodes; j++) {
            Sum += m_WSqrd(k,j) * m_pDotN(j);
        }
    }
    return Sum / m_sigmaQuad;
}

Eigen::VectorXd RBMJastrow::computeParameterGradient() {
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

void RBMJastrow::updateVectors() {
    m_v = m_b + m_W.transpose() * m_positions;
    Eigen::VectorXd m_e = m_v.array().exp();
    m_p = (m_e + Eigen::VectorXd::Ones(m_numberOfHiddenNodes)).cwiseInverse();
    m_n = m_e.cwiseProduct(m_p);
    m_pDotN   = m_p.cwiseProduct(m_n);
}

void RBMJastrow::updateRatio() {
    double Prod = 1;
    for(int j=0; j<m_numberOfHiddenNodes; j++) {
        Prod *= m_pOld(j)/m_p(j);
    }
    m_probabilityRatio  = Prod * Prod;
}
