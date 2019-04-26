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

double SimpleJastrow::calculateDistanceMatrixElement(const int i, const int j) {
    double dist = 0;
    int parti   = m_numberOfDimensions*i;
    int partj   = m_numberOfDimensions*j;
    for(int d=0; d<m_numberOfDimensions; d++) {
        double diff = m_positions(parti+d)-m_positions(partj+d);
        dist += diff*diff;
    }
    return sqrt(dist);
}

void SimpleJastrow::calculateDistanceMatrix() {
    m_distanceMatrix = Eigen::MatrixXd::Zero(m_numberOfParticles, m_numberOfParticles);
    for(int i=0; i<m_numberOfParticles; i++) {
        for(int j=i+1; j<m_numberOfParticles; j++) {
            m_distanceMatrix(i,j) = calculateDistanceMatrixElement(i,j);
            m_distanceMatrix(j,i) = m_distanceMatrix(i,j);
        }
    }
}

void SimpleJastrow::calculateDistanceMatrixCross(const int par) {
    for(int i=0; i<m_numberOfParticles; i++) {
        m_distanceMatrix(par, i) = calculateDistanceMatrixElement(par, i);
        m_distanceMatrix(i, par) = m_distanceMatrix(par, i);
    }
}

void SimpleJastrow::calculateG(int pRand) {
    for(int i=0; i<m_numberOfFreeDimensions; i++) {
        m_g(pRand,i) = m_positions(pRand) - m_positions(i);
        m_g(i,pRand) = -m_g(pRand,i);
    }
}

void SimpleJastrow::initializeArrays(Eigen::VectorXd positions) {
    m_positions = positions;
    m_probabilityRatio  = 1;
    calculateDistanceMatrix();
    m_g     = Eigen::MatrixXd::Zero(m_numberOfFreeDimensions, m_numberOfFreeDimensions);
    for(int i=0; i<m_numberOfFreeDimensions; i++) {
        for(int j=i; j<m_numberOfFreeDimensions; j++) {
            m_g(i,j) = m_positions(i) - m_positions(j);
            m_g(j,i) = -m_g(i,j);
        }
    }
    setArrays();
}

void SimpleJastrow::updateArrays(const Eigen::VectorXd positions, const int changedCoord) {
    int particle = int(changedCoord/m_numberOfDimensions);
    setArrays();

    m_positions     = positions;

    calculateDistanceMatrixCross(particle);
    calculateProbabilityRatio(particle);
    calculateG(changedCoord);
}

void SimpleJastrow::calculateProbabilityRatio(int particle) {
    double ratio = 0;
    for(int i=particle; i<m_numberOfParticles; i++) {
        ratio += m_beta(particle, i) * (m_distanceMatrix(particle, i) - m_distanceMatrixOld(particle, i));
    }
    m_probabilityRatio = exp(2*ratio);
}

void SimpleJastrow::setArrays() {
    m_positionsOld  = m_positions;
    m_distanceMatrixOld = m_distanceMatrix;
    m_probabilityRatioOld = m_probabilityRatio;
}

void SimpleJastrow::resetArrays() {
    m_positions         = m_positionsOld;
    m_distanceMatrix    = m_distanceMatrixOld;
    m_probabilityRatio  = m_probabilityRatioOld;
}

void SimpleJastrow::updateParameters(const Eigen::MatrixXd parameters, const int elementNumber) {
    m_elementNumber                     = elementNumber;
    m_maxNumberOfParametersPerElement   = m_system->getMaxNumberOfParametersPerElement();
    Eigen::VectorXd betaFlatten = parameters.row(m_elementNumber).head(m_numberOfFreeDimensions*m_numberOfFreeDimensions);
    Eigen::Map<Eigen::MatrixXd> beta(betaFlatten.data(), m_numberOfFreeDimensions, m_numberOfFreeDimensions);
    m_beta     = beta;
}

double SimpleJastrow::evaluateRatio() {
    return m_probabilityRatio;
}

double SimpleJastrow::computeGradient(const int k) {
    int k_p = int(k/m_numberOfDimensions);  //Particle associated with k
    int k_d = k%m_numberOfDimensions;       //Dimension associated with k

    double derivative = 0;
    for(int j_p=0; j_p<m_numberOfParticles; j_p++) {
        int j = j_p * m_numberOfDimensions + k_d;
        if(j_p!=k_p) {
            derivative += m_beta(k_p,j_p) * m_g(k,j) / m_distanceMatrix(k_p,j_p);
        }
    }
    return derivative;
}

double SimpleJastrow::computeLaplacian() {
    double derivative = 0;
    for(int i=0; i<m_numberOfFreeDimensions; i++) {
        int i_p = int(i/m_numberOfDimensions);  //Particle associated with k
        int i_d = i%m_numberOfDimensions;       //Dimension associated with k
        for(int j_p=0; j_p<m_numberOfParticles; j_p++) {
            int j = j_p * m_numberOfDimensions + i_d;
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
