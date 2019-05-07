#include "simplejastrow.h"
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include <iostream>

SimpleJastrow::SimpleJastrow(System* system) :
        WaveFunction(system) {
    m_numberOfParticles                 = m_system->getNumberOfParticles();
    m_numberOfDimensions                = m_system->getNumberOfDimensions();
    m_numberOfFreeDimensions            = m_system->getNumberOfFreeDimensions();
    m_numberOfParameters                = m_numberOfParticles * m_numberOfParticles;
}

void SimpleJastrow::calculateG(const unsigned int changedCoord) {
    for(unsigned int i=0; i<m_numberOfFreeDimensions; i++) {
        m_g(changedCoord,i) = m_positions(changedCoord) - m_positions(i);
        m_g(i,changedCoord) = -m_g(changedCoord,i);
    }
}

void SimpleJastrow::initializeArrays(const Eigen::VectorXd positions, const Eigen::VectorXd radialVector, const Eigen::MatrixXd distanceMatrix) {
    m_positions         = positions;
    m_distanceMatrix    = distanceMatrix;
    m_probabilityRatio  = 1;
    m_g     = Eigen::MatrixXd::Zero(m_numberOfFreeDimensions, m_numberOfFreeDimensions);
    for(unsigned int i=0; i<m_numberOfFreeDimensions; i++) {
        for(unsigned int j=i; j<m_numberOfFreeDimensions; j++) {
            m_g(i,j) = m_positions(i) - m_positions(j);
            m_g(j,i) = -m_g(i,j);
        }
    }
    setArrays();
}

void SimpleJastrow::updateArrays(const Eigen::VectorXd positions, const Eigen::VectorXd radialVector, const Eigen::MatrixXd distanceMatrix, const unsigned int changedCoord) {
    unsigned int particle = (unsigned int)(changedCoord/m_numberOfDimensions);
    setArrays();

    m_positions         = positions;
    m_distanceMatrix    = distanceMatrix;

    calculateProbabilityRatio(particle);
    calculateG(changedCoord);
}

void SimpleJastrow::calculateProbabilityRatio(unsigned int particle) {
    double ratio = 0;
    for(unsigned int i=particle; i<m_numberOfParticles; i++) {
        ratio += m_beta(particle, i) * (m_distanceMatrix(particle, i) - m_distanceMatrixOld(particle, i));
    }
    m_probabilityRatio = exp(2*ratio);
}

void SimpleJastrow::setArrays() {
    m_positionsOld          = m_positions;
    m_distanceMatrixOld     = m_distanceMatrix;
    m_probabilityRatioOld   = m_probabilityRatio;
}

void SimpleJastrow::resetArrays() {
    m_positions             = m_positionsOld;
    m_distanceMatrix        = m_distanceMatrixOld;
    m_probabilityRatio      = m_probabilityRatioOld;
}

void SimpleJastrow::updateParameters(const Eigen::MatrixXd parameters, const unsigned short elementNumber) {
    m_elementNumber                     = elementNumber;
    m_maxNumberOfParametersPerElement   = m_system->getMaxNumberOfParametersPerElement();
    Eigen::VectorXd betaFlatten = parameters.row(m_elementNumber).head(m_numberOfFreeDimensions*m_numberOfFreeDimensions);
    Eigen::Map<Eigen::MatrixXd> beta(betaFlatten.data(), m_numberOfFreeDimensions, m_numberOfFreeDimensions);
    m_beta     = beta;
}

double SimpleJastrow::evaluateRatio() {
    return m_probabilityRatio;
}

double SimpleJastrow::computeGradient(const unsigned int k) {
    unsigned int k_p = (unsigned int)(k/m_numberOfDimensions);  //Particle associated with k
    unsigned int k_d = k%m_numberOfDimensions;       //Dimension associated with k

    double derivative = 0;
    for(unsigned int j_p=0; j_p<m_numberOfParticles; j_p++) {
        unsigned int j = j_p * m_numberOfDimensions + k_d;
        if(j_p!=k_p) {
            derivative += m_beta(k_p,j_p) * m_g(k,j) / m_distanceMatrix(k_p,j_p);
        }
    }
    return derivative;
}

double SimpleJastrow::computeLaplacian() {
    double derivative = 0;
    for(unsigned int i=0; i<m_numberOfFreeDimensions; i++) {
        unsigned int i_p = (unsigned int)(i/m_numberOfDimensions);  //Particle associated with k
        unsigned int i_d = i%m_numberOfDimensions;       //Dimension associated with k
        for(unsigned int j_p=0; j_p<m_numberOfParticles; j_p++) {
            unsigned int j = j_p * m_numberOfDimensions + i_d;
            if(j_p!=i_p) {
                derivative += m_beta(i_p,j_p) * (1-m_g(i,j)*m_g(i,j)) / m_distanceMatrix(i_p,j_p);
            }
        }
    }
    return derivative;
}

Eigen::VectorXd SimpleJastrow::computeParameterGradient() {
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);
    Eigen::Map<Eigen::VectorXd> gradients2(m_distanceMatrix.data(), m_numberOfParticles*m_numberOfParticles);
    gradients.head(m_numberOfParticles*m_numberOfParticles) = gradients2;
    return gradients;
}
