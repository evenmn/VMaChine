#include "nqsjastrowold.h"
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include <iostream>

NQSJastrowOld::NQSJastrowOld(System* system) :
        WaveFunction(system) {
    m_numberOfHiddenNodes               = m_system->getNumberOfHiddenNodes();
    m_numberOfFreeDimensions            = m_system->getNumberOfFreeDimensions();
    m_maxNumberOfParametersPerElement   = m_system->getMaxNumberOfParametersPerElement();
    double sigma                        = m_system->getWidth();
    m_sigmaSqrd                         = sigma*sigma;
}

void NQSJastrowOld::updateArrays(const Eigen::VectorXd positions, const int pRand) {
    m_oldPositions = m_positions;
    m_positions = positions;

    m_oldV = m_v;
    m_oldN = m_n;
    m_oldP = m_p;

    m_v = m_b + m_W.transpose() * positions;
    for(int i=0; i<m_numberOfHiddenNodes; i++) {
        m_n(i) = 1/(1 + exp(-m_v(i)));
        m_p(i) = 1/(1 + exp(+m_v(i)));
    }

    m_oldRatio = m_ratio;
    double Prod = 1;
    for(int j=0; j<m_numberOfHiddenNodes; j++) {
        Prod *= (1 + exp(-m_v(j))) / (exp(-m_v(j)) + exp((m_oldPositions(pRand) - m_positions(pRand)) * m_W(pRand,j) / m_sigmaSqrd));
    }
    m_ratio = Prod * Prod;
}

void NQSJastrowOld::resetArrays(int pRand) {
    m_positions = m_oldPositions;
    m_v         = m_oldV;
    m_n         = m_oldN;
    m_p         = m_oldP;
    m_ratio     = m_oldRatio;
}

void NQSJastrowOld::initializeArrays(const Eigen::VectorXd positions) {
    m_positions = positions;
    m_v = m_b + m_W.transpose() * positions;

    m_n = Eigen::VectorXd::Zero(m_numberOfHiddenNodes);
    m_p = Eigen::VectorXd::Zero(m_numberOfHiddenNodes);
    for(int i=0; i<m_numberOfHiddenNodes; i++) {
        m_n(i) = 1/(1 + exp(-m_v(i)));
        m_p(i) = 1/(1 + exp(+m_v(i)));
    }

    m_ratio = 1;
}

void NQSJastrowOld::updateParameters(Eigen::MatrixXd parameters, const int elementNumber) {
    m_elementNumber = elementNumber;
    Eigen::VectorXd XXX = parameters.row(m_elementNumber).segment(m_numberOfHiddenNodes, m_numberOfFreeDimensions*m_numberOfHiddenNodes);
    Eigen::Map<Eigen::MatrixXd> W(XXX.data(), m_numberOfFreeDimensions, m_numberOfHiddenNodes);
    m_W = W;

    m_b = parameters.row(m_elementNumber).head(m_numberOfHiddenNodes);
}

int fromWToParameterIndex2(int i, int j, int numberOfFreeDimensions) {
    return j*numberOfFreeDimensions + i;
}

double NQSJastrowOld::evaluateRatio() {
    return m_ratio;
}

double NQSJastrowOld::computeFirstDerivative(const int k) {
    double Sum = 0;
    for(int j=0; j<m_numberOfHiddenNodes; j++) {
        Sum += m_W(k,j) * m_n(j);
    }
    return Sum / m_sigmaSqrd;
}

double NQSJastrowOld::computeSecondDerivative() {
    double Sum = 0;
    for(int k=0; k<m_numberOfFreeDimensions; k++) {
        for(int j=0; j<m_numberOfHiddenNodes; j++) {
            Sum += m_W(k,j) * m_W(k,j) * m_n(j) * m_p(j);
        }
    }
    return Sum / (m_sigmaSqrd * m_sigmaSqrd);
}

Eigen::VectorXd NQSJastrowOld::computeFirstEnergyDerivative(const int k) {
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);

    // Update b
    for(int j=0; j<m_numberOfHiddenNodes; j++) {
        gradients(j) = m_W(k,j) * m_n(j) * m_p(j) / m_sigmaSqrd;
    }


    // Update W
    for(int l=0; l<m_numberOfFreeDimensions; l++) {
        for(int m=0; m<m_numberOfHiddenNodes; m++) {
            int n = fromWToParameterIndex2(l, m, m_numberOfFreeDimensions);
            if(l == k) {
                gradients(n + m_numberOfHiddenNodes) = m_n(m) * (m_sigmaSqrd + m_W(k,m) * m_p(m) * m_positions(k));
            }
            else {
                gradients(n + m_numberOfHiddenNodes) = m_W(k, m) * m_n(m) * m_p(m) * m_positions(l);
            }
        }
    }
    return -0.5 * gradients / (m_sigmaSqrd * m_sigmaSqrd);
}

Eigen::VectorXd NQSJastrowOld::computeSecondEnergyDerivative() {
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);

    // Update b
    for(int k=0; k<m_numberOfFreeDimensions; k++) {
        for(int j=0; j<m_numberOfHiddenNodes; j++) {
            gradients(j) = m_W(k,j) * m_W(k,j) * m_n(j) * m_p(j) * (m_p(j) - m_n(j));
        }
    }

    // Update W
    for(int l=0; l<m_numberOfFreeDimensions; l++) {
        for(int m=0; m<m_numberOfHiddenNodes; m++) {
            int k = fromWToParameterIndex2(l, m, m_numberOfFreeDimensions);
            gradients(k + m_numberOfHiddenNodes) = m_n(m) * m_p(m) * (2 * m_W(l,m) + m_W.cwiseAbs2().colwise().sum()(m) * m_positions(l) * (m_p(m) - m_n(m)));
        }
    }

    return -0.5 * gradients / (m_sigmaSqrd * m_sigmaSqrd);
}
