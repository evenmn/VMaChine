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
}

void NQSJastrow::updateArrays(const Eigen::VectorXd positions, const int pRand) {
    m_positionsOld = m_positions;
    m_positions    = positions;

    m_vOld = m_v;
    m_nOld = m_n;
    m_pOld = m_p;

    m_v = m_b + m_W.transpose() * positions;
    for(int i=0; i<m_numberOfHiddenNodes; i++) {
        m_n(i) = 1/(1 + exp(-m_v(i)));
        m_p(i) = 1/(1 + exp(+m_v(i)));
    }

    m_probabilityRatioOld = m_probabilityRatio;
    double Prod = 1;
    for(int j=0; j<m_numberOfHiddenNodes; j++) {
        //Prod *= (1 + exp(-m_v(j))) / (exp(-m_v(j)) + exp((m_oldPositions(pRand) - m_positions(pRand)) * m_W(pRand,j) / m_sigmaSqrd));
        Prod *= (1 + exp(m_v(j)))/(1 + exp(m_vOld(j)));
    }
    m_probabilityRatio  = Prod * Prod;
}

void NQSJastrow::resetArrays(int pRand) {
    m_positions         = m_positionsOld;
    m_v                 = m_vOld;
    m_n                 = m_nOld;
    m_p                 = m_pOld;
    m_probabilityRatio  = m_probabilityRatioOld;
}

void NQSJastrow::initializeArrays(const Eigen::VectorXd positions) {
    m_positions = positions;
    m_v = m_b + m_W.transpose() * positions;

    m_n = Eigen::VectorXd::Zero(m_numberOfHiddenNodes);
    m_p = Eigen::VectorXd::Zero(m_numberOfHiddenNodes);
    for(int i=0; i<m_numberOfHiddenNodes; i++) {
        m_n(i) = 1/(1 + exp(-m_v(i)));
        m_p(i) = 1/(1 + exp(+m_v(i)));
    }

    m_probabilityRatio = 1;
}

void NQSJastrow::updateParameters(Eigen::MatrixXd parameters, const int elementNumber) {
    m_elementNumber = elementNumber;
    Eigen::VectorXd XXX = parameters.row(m_elementNumber).segment(m_numberOfHiddenNodes, m_numberOfFreeDimensions*m_numberOfHiddenNodes);
    Eigen::Map<Eigen::MatrixXd> W(XXX.data(), m_numberOfFreeDimensions, m_numberOfHiddenNodes);
    m_W = W;
    m_b = parameters.row(m_elementNumber).head(m_numberOfHiddenNodes);
}

int fromWToParameterIndex(int i, int j, int numberOfFreeDimensions) {
    return j*numberOfFreeDimensions + i;
}

double NQSJastrow::evaluateRatio() {
    return m_probabilityRatio;
}

double NQSJastrow::computeFirstDerivative(const int k) {
    double Sum = 0;
    for(int j=0; j<m_numberOfHiddenNodes; j++) {
        Sum += m_W(k,j) * m_n(j);
    }
    //return double(m_W.row(k) * m_n)/sigma_sqrd;
    return Sum / m_sigmaSqrd;
}

double NQSJastrow::computeSecondDerivative() {
    double Sum = 0;
    for(int k=0; k<m_numberOfFreeDimensions; k++) {
        for(int j=0; j<m_numberOfHiddenNodes; j++) {
            Sum += m_W(k,j) * m_W(k,j) * m_p(j) * m_n(j);
        }
    }
    //return (m_W.cwiseAbs2()*m_p.cwiseProduct(m_n)).sum();
    return Sum / (m_sigmaSqrd * m_sigmaSqrd);
}

Eigen::VectorXd NQSJastrow::computeFirstEnergyDerivative(const int k) {
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);

    // Update b
    for(int j=0; j<m_numberOfHiddenNodes; j++) {
        gradients(j) = m_W(k,j) * m_p(j) * m_n(j) / m_sigmaSqrd;
    }

    // Update W
    for(int l=0; l<m_numberOfFreeDimensions; l++) {
        for(int m=0; m<m_numberOfHiddenNodes; m++) {
            int n = fromWToParameterIndex(l, m, m_numberOfFreeDimensions);
            gradients(n + m_numberOfHiddenNodes) = m_W(k,m) * m_p(m) * m_n(m) * m_positions(l);
            if(l == k) {
                gradients(n + m_numberOfHiddenNodes) += m_n(m) * m_sigmaSqrd;
            }
        }
    }
    return -0.5 * gradients / (m_sigmaSqrd * m_sigmaSqrd);
}

Eigen::VectorXd NQSJastrow::computeSecondEnergyDerivative() {
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);

    // Update b
    for(int k=0; k<m_numberOfFreeDimensions; k++) {
        for(int j=0; j<m_numberOfHiddenNodes; j++) {
            gradients(j) = m_W(k,j) * m_W(k,j) * m_n(j) * m_p(j) * (m_p(j) - m_n(j));
        }
    }

    // Update W
    for(int k=0; k<m_numberOfFreeDimensions; k++) {
        for(int l=0; l<m_numberOfFreeDimensions; l++) {
            for(int m=0; m<m_numberOfHiddenNodes; m++) {
                int n = fromWToParameterIndex(l, m, m_numberOfFreeDimensions);
                gradients(n + m_numberOfHiddenNodes) = m_W(k,m) * m_W(k,m) * m_positions(l) * m_p(m) * m_n(m) * (m_p(m) - m_n(m)) / m_sigmaSqrd;
                if(l == k) {
                    gradients(n + m_numberOfHiddenNodes) += 2 * m_W(k,m) * m_p(m) * m_n(m);
                }
            }
        }
    }

    return -0.5 * gradients / (m_sigmaSqrd * m_sigmaSqrd);
}
