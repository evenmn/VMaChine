#include "samsethjastrow.h"
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include <iostream>

SamsethJastrow::SamsethJastrow(System* system) :
        WaveFunction(system) {
    m_numberOfParticles                 = m_system->getNumberOfParticles();
    m_numberOfDimensions                = m_system->getNumberOfDimensions();
    m_numberOfFreeDimensions            = m_system->getNumberOfFreeDimensions();
}

void SamsethJastrow::calculateG(int pRand) {
    for(int i=0; i<m_numberOfFreeDimensions; i++) {
        m_g(pRand,i) = m_positions(pRand) - m_positions(i);
        m_g(i,pRand) = -m_g(pRand,i);
    }
    m_gSqrd = m_g.cwiseAbs2();
}

void SamsethJastrow::calculateProbabilityRatio(int particle) {
    //double ratio = (m_beta.row(particle).transpose() * (m_h.row(particle) - m_hOld.row(particle))).sum();

    double ratio = 0;
    for(int i=particle; i<m_numberOfParticles; i++) {
        ratio += - m_betaSqrd * m_distanceMatrixSqrd(particle,i)/2    + m_betagamma*m_distanceMatrix(particle,i)  -  \
                 - m_betaSqrd * m_distanceMatrixSqrdOld(particle,i)/2 + m_betagamma*m_distanceMatrixOld(particle,i);
    }
    m_probabilityRatio = exp(ratio);
}

void SamsethJastrow::initializeArrays(const Eigen::VectorXd positions, const Eigen::VectorXd radialVector, const Eigen::MatrixXd distanceMatrix) {
    m_positions         = positions;
    m_distanceMatrix    = distanceMatrix;
    m_probabilityRatio  = 1;
    m_g     = Eigen::MatrixXd::Zero(m_numberOfFreeDimensions, m_numberOfFreeDimensions);
    for(int i=0; i<m_numberOfFreeDimensions; i++) {
        for(int j=i; j<m_numberOfFreeDimensions; j++) {
            m_g(i,j) = m_positions(i) - m_positions(j);
            m_g(j,i) = -m_g(i,j);
        }
    }
    setArrays();
}

void SamsethJastrow::updateArrays(const Eigen::VectorXd positions, const Eigen::VectorXd radialVector, const Eigen::MatrixXd distanceMatrix, const int changedCoord) {
    int particle = int(changedCoord/m_numberOfDimensions);
    setArrays();
    m_positions             = positions;
    m_distanceMatrix        = distanceMatrix;
    calculateG                  (changedCoord);
    calculateProbabilityRatio   (particle);
}

void SamsethJastrow::setArrays() {
    m_positionsOld          = m_positions;
    m_distanceMatrixOld     = m_distanceMatrix;
    m_distanceMatrixSqrdOld = m_distanceMatrixSqrd;
    m_gOld                  = m_g;
    m_gSqrdOld              = m_gSqrd;
    m_probabilityRatioOld   = m_probabilityRatio;
}

void SamsethJastrow::resetArrays() {
    m_positions             = m_positionsOld;
    m_distanceMatrix        = m_distanceMatrixOld;
    m_distanceMatrixSqrd    = m_distanceMatrixSqrdOld;
    m_g                     = m_gOld;
    m_gSqrd                 = m_gSqrdOld;
    m_probabilityRatio      = m_probabilityRatioOld;
}

void SamsethJastrow::updateParameters(const Eigen::MatrixXd parameters, const int elementNumber) {
    m_elementNumber                     = elementNumber;
    m_maxNumberOfParametersPerElement   = m_system->getMaxNumberOfParametersPerElement();
    m_beta                              = parameters(m_elementNumber, 0);
    m_gamma                             = parameters(m_elementNumber, 1);
    m_betagamma                         = fabs(m_beta * m_gamma);
    m_betaSqrd                          = m_beta * m_beta;
}

double SamsethJastrow::evaluateRatio() {
    return m_probabilityRatio;
}

double SamsethJastrow::computeGradient(const int k) {
    int k_p = int(k/m_numberOfDimensions);  //Particle associated with k
    int k_d = k%m_numberOfDimensions;       //Dimension associated with k

    double derivative = 0;
    for(int j_p=0; j_p<m_numberOfParticles; j_p++) {
        int j = j_p * m_numberOfDimensions + k_d;
        if(j_p!=k_p) {
            derivative += (m_betagamma / m_distanceMatrix(k_p,j_p) - m_betaSqrd) * m_g(k,j);
        }
    }
    return 2*derivative;
}

double SamsethJastrow::computeLaplacian() {
    double derivative = 0;
    for(int i=0; i<m_numberOfFreeDimensions; i++) {
        int i_p = int(i/m_numberOfDimensions);  //Particle associated with k
        int i_d = i%m_numberOfDimensions;       //Dimension associated with k
        for(int j_p=i_p+1; j_p<m_numberOfParticles; j_p++) {
            int j = j_p * m_numberOfDimensions + i_d;
            derivative += -m_betaSqrd * m_gSqrd(i,j) / m_distanceMatrixSqrd(i_p,j_p) + (m_betagamma/m_distanceMatrix(i_p,j_p) - m_betaSqrd) * (1 + m_gSqrd(i,j) / m_distanceMatrixSqrd(i_p,j_p));
        }
    }
    return 2 * derivative;
}

Eigen::VectorXd SamsethJastrow::computeParameterGradient() {
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);

    double updateBeta  = 0;
    double updateGamma = 0;
    for(int k=0; k<m_numberOfFreeDimensions; k++) {
        int k_p = int(k/m_numberOfDimensions);  //Particle associated with k
        for(int j_p=k_p+1; j_p<m_numberOfParticles; j_p++) {
            updateBeta  += -m_beta * m_distanceMatrixSqrd(k_p, j_p) + m_gamma * m_distanceMatrix(k_p,j_p);
            updateGamma += m_beta * m_distanceMatrix(k_p, j_p);
        }
    }
    gradients(0) = updateBeta;
    gradients(1) = updateGamma;
    return gradients;
}
