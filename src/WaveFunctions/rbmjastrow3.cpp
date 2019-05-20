#include "rbmjastrow3.h"
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include <iostream>

RBMJastrow3::RBMJastrow3(System* system) :
        WaveFunction(system) {
    m_numberOfParticles                 = m_system->getNumberOfParticles();
    m_numberOfDimensions                = m_system->getNumberOfDimensions();
    m_numberOfHiddenNodes               = m_system->getNumberOfHiddenNodes();
    m_numberOfFreeDimensions            = m_system->getNumberOfFreeDimensions();
    m_numberOfParameters                = m_numberOfParticles * m_numberOfParticles * m_numberOfHiddenNodes + m_numberOfHiddenNodes;
    double sigma                        = 1; //m_system->getWidth();
    m_sigmaSqrd                         = sigma*sigma;
    m_sigmaQuad                         = m_sigmaSqrd*m_sigmaSqrd;
}

void RBMJastrow3::setConstants(const unsigned int elementNumber) {
    m_maxNumberOfParametersPerElement   = m_system->getMaxParameters();
    m_elementNumber                     = elementNumber;
}

void RBMJastrow3::calculateG(int changedCoord) {
    for(int i=changedCoord+1; i<m_numberOfFreeDimensions; i++) {
        m_g(changedCoord,i) = m_positions(changedCoord) - m_positions(i);
    }
    for(int i=0; i<changedCoord; i++) {
        m_g(i,changedCoord) = m_positions(i) - m_positions(changedCoord);
    }
}

void RBMJastrow3::updateVectors() {
    for(int l=0; l<m_numberOfHiddenNodes; l++) {
        m_v(l) = m_b(l);
        for(int i=0; i<m_numberOfParticles; i++) {
            for(int j=i+1; j<m_numberOfParticles; j++) {
                m_v(l) += m_W(l*i,j) * m_distanceMatrix(i,j) / m_sigmaSqrd;
            }
        }
    }

    Eigen::VectorXd m_e = m_v.array().exp();
    m_p = (m_e + Eigen::VectorXd::Ones(m_numberOfHiddenNodes)).cwiseInverse();
    m_n = m_e.cwiseProduct(m_p);
}

void RBMJastrow3::updateRatio() {
    double Prod = 1;
    for(int j=0; j<m_numberOfHiddenNodes; j++) {
        Prod *= m_pOld(j)/m_p(j);
    }
    m_probabilityRatio  = Prod * Prod;
}

void RBMJastrow3::updateArrays(const Eigen::VectorXd positions, const Eigen::VectorXd radialVector, const Eigen::MatrixXd distanceMatrix, const int changedCoord) {
    int particle = int(changedCoord/m_numberOfDimensions);
    m_positions = positions;
    m_distanceMatrix = distanceMatrix;
    calculateG(changedCoord);
    updateVectors();
    updateRatio();
}

void RBMJastrow3::setArrays() {
    m_positionsOld          = m_positions;
    m_vOld                  = m_v;
    m_nOld                  = m_n;
    m_pOld                  = m_p;
    m_gOld                  = m_g;
    m_distanceMatrixOld     = m_distanceMatrix;
    m_probabilityRatioOld   = m_probabilityRatio;
}

void RBMJastrow3::resetArrays() {
    m_positions             = m_positionsOld;
    m_v                     = m_vOld;
    m_n                     = m_nOld;
    m_p                     = m_pOld;
    m_g                     = m_gOld;
    m_distanceMatrix        = m_distanceMatrixOld;
    m_probabilityRatio      = m_probabilityRatioOld;
}

void RBMJastrow3::initializeArrays(const Eigen::VectorXd positions, const Eigen::VectorXd radialVector, const Eigen::MatrixXd distanceMatrix) {
    m_positions         = positions;
    m_distanceMatrix    = distanceMatrix;
    m_probabilityRatio  = 1;

    m_v  = Eigen::VectorXd::Zero(m_numberOfHiddenNodes);
    m_n  = Eigen::VectorXd::Zero(m_numberOfHiddenNodes);
    m_p  = Eigen::VectorXd::Zero(m_numberOfHiddenNodes);
    m_g  = Eigen::MatrixXd::Zero(m_numberOfFreeDimensions, m_numberOfFreeDimensions);
    for(int i=0; i<m_numberOfFreeDimensions; i++) {
        for(int j=i; j<m_numberOfFreeDimensions; j++) {
            m_g(i,j) = m_positions(i) - m_positions(j);
            m_g(j,i) = -m_g(i,j);
        }
    }
    updateVectors();
}

void RBMJastrow3::updateParameters(Eigen::MatrixXd parameters) {
    Eigen::VectorXd wFlatten = parameters.row(m_elementNumber).segment(m_numberOfHiddenNodes, m_numberOfParticles*m_numberOfParticles*m_numberOfHiddenNodes);
    Eigen::Map<Eigen::MatrixXd> W(wFlatten.data(), m_numberOfParticles*m_numberOfHiddenNodes, m_numberOfParticles);
    m_W     = W;
    m_WSqrd = W.cwiseAbs2();
    m_b     = parameters.row(m_elementNumber).head(m_numberOfHiddenNodes);
}

double RBMJastrow3::evaluateRatio() {
    return m_probabilityRatio;
}

double RBMJastrow3::computeGradient(const int k) {
    int k_p = int(k/m_numberOfDimensions);  //Particle associated with k
    int k_d = k%m_numberOfDimensions;       //Dimension associated with k
    double derivative = 0;
    for(int l=0; l<m_numberOfHiddenNodes; l++) {
        double firstder = 0;
        for(int i_p=0; i_p<k_p; i_p++) {
            int i = i_p * m_numberOfDimensions + k_d;
            firstder -= m_W(l*i_p,k_p) * m_g(i,k) / m_distanceMatrix(i_p,k_p);
        }
        for(int i_p=k_p+1; i_p<m_numberOfParticles; i_p++) {
            int i = i_p * m_numberOfDimensions + k_d;
            firstder += m_W(l*k_p,i_p) * m_g(k,i) / m_distanceMatrix(k_p,i_p);
        }
        derivative += m_n(l) * firstder;
    }
    return derivative / m_sigmaSqrd;
}

double RBMJastrow3::computeLaplacian() {
    double derivative = 0;
    for(int k=0; k<m_numberOfFreeDimensions; k++) {
        int k_p = int(k/m_numberOfDimensions);  //Particle associated with k
        int k_d = k%m_numberOfDimensions;       //Dimension associated with k
        for(int l=0; l<m_numberOfHiddenNodes; l++) {
            double firstder = 0;
            double secondder = 0;
            for(int i_p=0; i_p<k_p; i_p++) {
                int i = i_p * m_numberOfDimensions + k_d;
                firstder -= m_W(l*i_p,k_p) * m_g(i,k) / m_distanceMatrix(i_p,k_p);
                secondder -= m_W(l*i_p,k_p) * (1 - m_g(i,k) * m_g(i,k) / (m_distanceMatrix(i_p,k_p) * m_distanceMatrix(i_p,k_p)));
            }
            for(int i_p=k_p+1; i_p<m_numberOfParticles; i_p++) {
                int i = i_p * m_numberOfDimensions + k_d;
                firstder += m_W(l*k_p,i_p) * m_g(k,i) / m_distanceMatrix(k_p,i_p);
                secondder += m_W(l*k_p,i_p) * (1 - m_g(k,i) * m_g(k,i) / (m_distanceMatrix(k_p,i_p) * m_distanceMatrix(k_p,i_p)));
            }
            derivative += m_n(l) * (secondder / m_sigmaSqrd + m_p(l) * firstder * firstder);
        }
    }
    return derivative / m_sigmaSqrd;
}

Eigen::VectorXd RBMJastrow3::computeParameterGradient() {
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);
    for(int l=0; l<m_numberOfHiddenNodes; l++) {
        gradients(l) = m_n(l);
    }
    for(int m=0; m<m_numberOfParticles; m++) {
        for(int n=m+1; n<m_numberOfParticles; n++) {
            for(int o=0; o<m_numberOfHiddenNodes; o++) {
                int p = o * m_numberOfParticles*m_numberOfHiddenNodes + m*m_numberOfParticles + n + m_numberOfHiddenNodes;
                gradients(p) = m_n(o) * m_distanceMatrix(m,n) / m_sigmaSqrd;
            }
        }
    }
    return gradients;
}
