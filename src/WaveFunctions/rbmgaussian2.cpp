#include "rbmgaussian2.h"
#include <cassert>
#include <iostream>
#include "../system.h"

RBMGaussian2::RBMGaussian2(System* system) :
        WaveFunction(system) {
    m_numberOfParticles                 = m_system->getNumberOfParticles();
    m_numberOfDimensions                = m_system->getNumberOfDimensions();
    m_numberOfFreeDimensions            = m_system->getNumberOfFreeDimensions();
    m_numberOfParameters                = m_numberOfParticles * m_numberOfParticles;
    m_omega                             = m_system->getFrequency();
    double sigma                        = m_system->getWidth();
    m_sigmaSqrd = sigma*sigma;
}

double RBMGaussian2::calculateDistanceMatrixElement(const int i, const int j) {
    double dist = 0;
    int parti   = m_numberOfDimensions*i;
    int partj   = m_numberOfDimensions*j;
    for(int d=0; d<m_numberOfDimensions; d++) {
        double diff = m_positions(parti+d)-m_positions(partj+d);
        dist += diff*diff;
    }
    return sqrt(dist);
}

void RBMGaussian2::calculateDistanceMatrix() {
    m_distanceMatrix = Eigen::MatrixXd::Zero(m_numberOfParticles, m_numberOfParticles);
    for(int i=0; i<m_numberOfParticles; i++) {
        for(int j=i+1; j<m_numberOfParticles; j++) {
            m_distanceMatrix(i,j) = calculateDistanceMatrixElement(i,j);
        }
    }
}

void RBMGaussian2::calculateDistanceMatrixCross(const int par) {
    for(int i=0; i<par; i++) {
        m_distanceMatrix(i, par) = calculateDistanceMatrixElement(i, par);
    }
    for(int j=par+1; j<m_numberOfParticles; j++) {
        m_distanceMatrix(par, j) = calculateDistanceMatrixElement(par, j);
    }
}

void RBMGaussian2::calculateG(int changedCoord) {
    for(int i=changedCoord+1; i<m_numberOfFreeDimensions; i++) {
        m_g(changedCoord,i) = m_positions(changedCoord) - m_positions(i);
    }
    for(int i=0; i<changedCoord; i++) {
        m_g(i,changedCoord) = m_positions(i) - m_positions(changedCoord);
    }
    m_gSqrd = m_g.cwiseAbs2();
}

void RBMGaussian2::updateParameters(const Eigen::MatrixXd parameters, const int elementNumber) {
    m_elementNumber                     = elementNumber;
    m_maxNumberOfParametersPerElement   = m_system->getMaxNumberOfParametersPerElement();
    Eigen::VectorXd aFlatten = parameters.row(m_elementNumber).head(m_numberOfParticles*m_numberOfParticles);
    Eigen::Map<Eigen::MatrixXd> a(aFlatten.data(), m_numberOfParticles, m_numberOfParticles);
    m_a     = a;
}

void RBMGaussian2::initializeArrays(const Eigen::VectorXd positions) {
    m_positions             = positions;

    calculateDistanceMatrix();

    m_Xa                    = m_distanceMatrix - m_a;
    m_probabilityRatio      = 1;

    m_g  = Eigen::MatrixXd::Zero(m_numberOfFreeDimensions, m_numberOfFreeDimensions);
    for(int i=0; i<m_numberOfFreeDimensions; i++) {
        for(int j=i+1; j<m_numberOfFreeDimensions; j++) {
            m_g(i,j) = m_positions(i) - m_positions(j);
        }
    }
    m_gSqrd = m_g.cwiseAbs2();

    setArrays();
}

void RBMGaussian2::updateArrays(const Eigen::VectorXd positions, const int changedCoord) {
    int particle = int(changedCoord/m_numberOfDimensions);
    setArrays();

    m_positions             = positions;
    calculateDistanceMatrixCross(particle);
    m_Xa                    = m_distanceMatrix - m_a;

    calculateG(changedCoord);

    double ratio = 0;
    for(int j=particle+1; j<m_numberOfParticles; j++) {
        ratio += m_XaOld(particle,j) - m_Xa(particle,j);
    }
    m_probabilityRatio = exp(ratio/m_sigmaSqrd);
}

void RBMGaussian2::setArrays() {
    m_positionsOld          = m_positions;
    m_distanceMatrixOld     = m_distanceMatrix;
    m_XaOld                 = m_Xa;
    m_gOld                  = m_g;
    m_gSqrdOld              = m_gSqrd;
    m_probabilityRatioOld   = m_probabilityRatio;
}

void RBMGaussian2::resetArrays() {
    m_positions             = m_positionsOld;
    m_distanceMatrix        = m_distanceMatrixOld;
    m_Xa                    = m_XaOld;
    m_g                     = m_gOld;
    m_gSqrd                 = m_gSqrdOld;
    m_probabilityRatio      = m_probabilityRatioOld;
}

double RBMGaussian2::evaluateRatio() {
    return m_probabilityRatio;
}

double RBMGaussian2::computeGradient(const int k) {
    int k_p = int(k/m_numberOfDimensions);  //Particle associated with k
    int k_d = k%m_numberOfDimensions;       //Dimension associated with k

    double sum = 0;
    for(int i_p=0; i_p<k_p; i_p++) {
        int i = i_p * m_numberOfDimensions + k_d;
        sum += m_Xa(i_p,k_p) * m_g(i,k);
    }
    for(int j_p = k_p+1; j_p<m_numberOfParticles; j_p++) {
        int j = j_p * m_numberOfDimensions + k_d;
        sum -= m_Xa(k_p,j_p) * m_g(k,j);
    }
    return sum/(2*m_sigmaSqrd);
}

double RBMGaussian2::computeLaplacian() {
    double sum = 0;
    for(int k=0; k<m_numberOfFreeDimensions; k++) {
        int k_p = int(k/m_numberOfDimensions);  //Particle associated with k
        int k_d = k%m_numberOfDimensions;       //Dimension associated with k

        for(int i_p=0; i_p<k_p; i_p++) {
         int i = i_p * m_numberOfDimensions + k_d;
         sum -= m_g(i,k);
         sum += m_Xa(i_p,k_p) * (1 - m_gSqrd(i,k)) / m_distanceMatrix(i_p,k_p);
        }
        for(int j_p = k_p+1; j_p<m_numberOfParticles; j_p++) {
         int j = j_p * m_numberOfDimensions + k_d;
         sum += m_g(k,j);
         sum -= m_Xa(k_p,j_p) * (1 - m_gSqrd(k,j)) / m_distanceMatrix(k_p,j_p);
        }
    }
    return -sum/(2*m_sigmaSqrd);
}

Eigen::VectorXd RBMGaussian2::computeParameterGradient() {
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);
    for(int i=0; i<m_numberOfParticles; i++) {
        for(int j=i+1; j<m_numberOfParticles; j++) {
            gradients(j*m_numberOfParticles+1) = m_Xa(i,j)/m_sigmaSqrd;
        }
    }
    return gradients;
}
