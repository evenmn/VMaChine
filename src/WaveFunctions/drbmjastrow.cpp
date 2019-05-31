#include "drbmjastrow.h"
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include <iostream>

DRBMJastrow::DRBMJastrow(System* system, int numberOfLayers) :
        WaveFunction(system) {
    m_numberOfLayers                    = numberOfLayers;
    m_numberOfHiddenNodes               = m_system->getNumberOfHiddenNodes();
    m_numberOfFreeDimensions            = m_system->getNumberOfFreeDimensions();
    m_numberOfParameters                = m_numberOfHiddenNodes*(1 + numberOfLayers*m_numberOfFreeDimensions);
    double sigma                        = 1; //m_system->getWidth();
    m_sigmaSqrd                         = sigma*sigma;
}

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

void DRBMJastrow::setConstants(const int elementNumber) {
    m_maxNumberOfParametersPerElement   = m_system->getMaxParameters();
    m_elementNumber                     = elementNumber;
}

void DRBMJastrow::updateGradient() {
    for(int k=0; k<m_numberOfFreeDimensions; k++) {
        for(int j=0; j<m_numberOfHiddenNodes; j++) {
            m_gradientPart(k,j) = 0;
            for(int n=1; n<m_numberOfLayers+1; n++) {
                m_gradientPart(k,j) += n*m_positionsPow(n,k)*m_W((n-1)*m_numberOfFreeDimensions+k,j) / pow(m_sigmaSqrd, 2*n);
            }
        }
    }
}

void DRBMJastrow::updateLaplacian() {
    for(int k=0; k<m_numberOfFreeDimensions; k++) {
        for(int j=0; j<m_numberOfHiddenNodes; j++) {
            m_laplacianPart(k,j) = 0;
            for(int n=0; n<m_numberOfLayers; n++) {
                m_gradientPart(k,j) += n*(n+1)*m_positionsPow(n,k)*m_W(n*m_numberOfFreeDimensions+k,j) / pow(m_sigmaSqrd, 2*(n+1));
            }
        }
    }
}

void DRBMJastrow::updateVectors() {
    m_v = m_b;
    for(int n=0; n<m_numberOfLayers; n++) {
        m_v += m_positionsPow.row(n+2)*m_W.block(n*m_numberOfFreeDimensions, 0, m_numberOfFreeDimensions, m_numberOfHiddenNodes) / pow(m_sigmaSqrd, 2*(n+1));
    }
    Eigen::VectorXd m_e = m_v.array().exp();
    m_p = (m_e + Eigen::VectorXd::Ones(m_numberOfHiddenNodes)).cwiseInverse();
    m_n = m_e.cwiseProduct(m_p);
}

void DRBMJastrow::updateRatio() {
    double Prod = 1;
    for(int j=0; j<m_numberOfHiddenNodes; j++) {
        Prod *= m_pOld(j)/m_p(j);
    }
    m_probabilityRatio  = Prod * Prod;
}

void DRBMJastrow::updateArrays(const Eigen::VectorXd positions, const Eigen::VectorXd radialVector, const Eigen::MatrixXd distanceMatrix, const int changedCoord) {
    m_positions     = positions;
    for(int n=1; n<m_numberOfLayers+1; n++) {
        m_positionsPow(n+1, changedCoord) = pow(m_positions(changedCoord),n);
    }
    updateVectors();
    updateRatio();
    updateGradient();
    updateLaplacian();
}

void DRBMJastrow::setArrays() {
    m_positionsOld          = m_positions;
    m_positionsPowOld       = m_positionsPow;
    m_gradientPartOld       = m_gradientPart;
    m_laplacianPartOld      = m_laplacianPart;
    m_vOld                  = m_v;
    m_nOld                  = m_n;
    m_pOld                  = m_p;
    m_probabilityRatioOld   = m_probabilityRatio;
}

void DRBMJastrow::resetArrays() {
    m_positions             = m_positionsOld;
    m_positionsPow          = m_positionsPowOld;
    m_gradientPart          = m_gradientPartOld;
    m_laplacianPart         = m_laplacianPartOld;
    m_v                     = m_vOld;
    m_n                     = m_nOld;
    m_p                     = m_pOld;
    m_probabilityRatio      = m_probabilityRatioOld;
}

void DRBMJastrow::initializeArrays(const Eigen::VectorXd positions, const Eigen::VectorXd radialVector, const Eigen::MatrixXd distanceMatrix) {
    m_positions         = positions;
    m_positionsPow      = Eigen::MatrixXd::Zero(m_numberOfLayers+2, m_numberOfFreeDimensions);
    for(int n=0; n<m_numberOfLayers+1; n++) {
        m_positionsPow.row(n+1) = m_positions.array().pow(n);
    }
    m_probabilityRatio  = 1;

    m_n             = Eigen::VectorXd::Zero(m_numberOfHiddenNodes);
    m_p             = Eigen::VectorXd::Zero(m_numberOfHiddenNodes);
    m_gradientPart  = Eigen::MatrixXd::Zero(m_numberOfFreeDimensions, m_numberOfHiddenNodes);
    m_laplacianPart = Eigen::MatrixXd::Zero(m_numberOfFreeDimensions, m_numberOfHiddenNodes);

    updateVectors();
    updateGradient();
    updateLaplacian();
}

void DRBMJastrow::updateParameters(Eigen::MatrixXd parameters) {
    m_b = parameters.row(m_elementNumber).head(m_numberOfHiddenNodes);
    m_W = Eigen::MatrixXd::Zero(m_numberOfLayers*m_numberOfFreeDimensions, m_numberOfHiddenNodes);
    for(int n=0; n<m_numberOfLayers; n++) {
        Eigen::VectorXd wFlatten = parameters.row(m_elementNumber).segment(m_numberOfHiddenNodes*(1+n*m_numberOfFreeDimensions), m_numberOfFreeDimensions*m_numberOfHiddenNodes);
        Eigen::Map<Eigen::MatrixXd> W(wFlatten.data(), m_numberOfFreeDimensions, m_numberOfHiddenNodes);
        m_W.block(n*m_numberOfFreeDimensions, 0, m_numberOfFreeDimensions, m_numberOfHiddenNodes) = W;
        //m_W.block(n*m_numberOfFreeDimensions, 0, m_numberOfFreeDimensions, m_numberOfHiddenNodes) = WaveFunction::reshape(wFlatten, m_numberOfFreeDimensions, m_numberOfHiddenNodes);
    }
    //m_W.block(m_numberOfFreeDimensions, 0, m_numberOfFreeDimensions, m_numberOfHiddenNodes) = Eigen::MatrixXd::Zero(m_numberOfFreeDimensions, m_numberOfHiddenNodes);
    //m_W.block(2*m_numberOfFreeDimensions, 0, m_numberOfFreeDimensions, m_numberOfHiddenNodes) = Eigen::MatrixXd::Zero(m_numberOfFreeDimensions, m_numberOfHiddenNodes);
    //std::cout << m_W << std::endl;
}

double DRBMJastrow::evaluateRatio() {
    return m_probabilityRatio;
}

double DRBMJastrow::computeGradient(const int k) {
    return m_gradientPart.row(k) * m_n;
}

double DRBMJastrow::computeLaplacian() {
    double sum = 0;
    for(int k=0; k<m_numberOfFreeDimensions; k++) {
        for(int j=0; j<m_numberOfHiddenNodes; j++) {
            sum += m_n(j) * (m_laplacianPart(k,j) + m_p(j) * m_gradientPart(k,j) * m_gradientPart(k,j));
        }
    }
    return sum;
}

Eigen::VectorXd DRBMJastrow::computeParameterGradient() {
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);

    for(int l=0; l<m_numberOfHiddenNodes; l++) {
        gradients(l) = m_n(l);
        for(int m=0; m<m_numberOfFreeDimensions; m++) {
            for(int n=0; n<m_numberOfLayers; n++) {
                int o = l * m_numberOfFreeDimensions + m + m_numberOfHiddenNodes * (1 + n*m_numberOfFreeDimensions);
                gradients(o) = m_positionsPow(n+2,m) * m_n(l) / pow(m_sigmaSqrd, 2*(n+1));
            }
        }
    }

    /*
    gradients.head(m_numberOfHiddenNodes) = m_n;
    for(int n=0; n<m_numberOfLayers; n++) {
        Eigen::MatrixXd out = m_positionsPow.row(n+2) * m_n.transpose();
        gradients.segment(m_numberOfHiddenNodes, m_numberOfHiddenNodes*m_numberOfFreeDimensions) = WaveFunction::flatten(out) / pow(m_sigmaSqrd, 2*(n+1));
    }
    */
    return gradients;
}
