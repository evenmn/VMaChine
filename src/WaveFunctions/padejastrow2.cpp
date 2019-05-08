#include "padejastrow2.h"
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include <iostream>

PadeJastrow2::PadeJastrow2(System* system) :
        WaveFunction(system) {
    m_numberOfParticles                 = m_system->getNumberOfParticles();
    m_numberOfDimensions                = m_system->getNumberOfDimensions();
    m_numberOfFreeDimensions            = m_system->getNumberOfFreeDimensions();
}

void PadeJastrow2::calculateF(const unsigned int particle) {
    for(unsigned int i=0; i<particle; i++) {
        m_f(i, particle) = 1/(1 + m_gamma * m_distanceMatrix(i, particle));
    }
    for(unsigned int j=particle+1; j<m_numberOfParticles; j++) {
        m_f(particle, j) = 1/(1 + m_gamma * m_distanceMatrix(particle, j));
    }
    m_fSqrd = m_f.cwiseAbs2();
}

void PadeJastrow2::calculateG(const unsigned int particle, const unsigned int changedCoord) {
    for(unsigned int i=0; i<changedCoord; i++) {
        unsigned int i_p = unsigned(i/m_numberOfDimensions);
        m_g(i, changedCoord) = (m_positions(i) - m_positions(changedCoord))/m_distanceMatrix(i_p, particle);
    }
    for(unsigned int j=changedCoord+1; j<m_numberOfFreeDimensions; j++) {
        unsigned int j_p = unsigned(j/m_numberOfDimensions);
        m_g(changedCoord, j) = (m_positions(changedCoord) - m_positions(j))/m_distanceMatrix(particle, j_p);
    }
    m_gSqrd = m_g.cwiseAbs2();
}

void PadeJastrow2::calculateH(const unsigned int particle) {
    for(unsigned int i=0; i<particle; i++) {
        m_h(i, particle) = m_distanceMatrix(i, particle) * m_f(i, particle);
    }
    for(unsigned int j=particle+1; j<m_numberOfParticles; j++) {
        m_h(particle, j) = m_distanceMatrix(particle, j) * m_f(particle, j);
    }
}

void PadeJastrow2::calculateProbabilityRatio(const unsigned int particle) {
    double ratio = 0;
    for(unsigned int j=particle+1; j<m_numberOfParticles; j++) {
        ratio += m_beta(particle, j) * (m_h(particle, j) - m_hOld(particle, j));
    }
    m_probabilityRatio = exp(2*ratio);
}


void PadeJastrow2::initializeBeta() {
    double symmetric, antisymmetric;
    if (m_numberOfDimensions == 2) {
        symmetric     = 1./3;
        antisymmetric = 1./1;
    }
    else if(m_numberOfDimensions == 3) {
        symmetric     = 1./4;
        antisymmetric = 1./2;
    }
    else if(m_numberOfDimensions == 4) {
        symmetric     = 1./5;
        antisymmetric = 1./3;
    }
    else {
        std::cout << "Number of dimensions should be in the range [1, 4]" << std::endl;
        exit(0);
    }

    m_beta = Eigen::MatrixXd::Zero(m_numberOfParticles, m_numberOfParticles);

    unsigned int numberOfParticlesHalf = m_numberOfParticles/2;
    for (unsigned int i = 0; i < m_numberOfParticles; i++) {
        for (unsigned int j = 0; j < m_numberOfParticles; j++) {

            if ((j <  numberOfParticlesHalf && i <  numberOfParticlesHalf) ||  \
                (j >= numberOfParticlesHalf && i >= numberOfParticlesHalf)) {
                m_beta(i, j) = symmetric;
            } else {
                m_beta(i, j) = antisymmetric;
            }
        }
    }
}

void PadeJastrow2::initializeArrays(const Eigen::VectorXd positions, const Eigen::VectorXd radialVector, const Eigen::MatrixXd distanceMatrix) {
    m_positions         = positions;
    m_distanceMatrix    = distanceMatrix;
    m_probabilityRatio  = 1;
    m_f     = (Eigen::MatrixXd::Ones(m_numberOfParticles, m_numberOfParticles) + m_gamma * m_distanceMatrix).cwiseInverse();
    m_g     = Eigen::MatrixXd::Zero(m_numberOfFreeDimensions, m_numberOfFreeDimensions);
    for(unsigned int i=0; i<m_numberOfFreeDimensions; i++) {
        for(unsigned int j=i+1; j<m_numberOfFreeDimensions; j++) {
            m_g(i,j) = m_positions(i) - m_positions(j);
        }
    }
    m_h = m_distanceMatrix.cwiseProduct(m_f);
    m_fSqrd = m_f.cwiseAbs2();
    m_gSqrd = m_g.cwiseAbs2();

    initializeBeta();
}

void PadeJastrow2::updateArrays(const Eigen::VectorXd positions, const Eigen::VectorXd radialVector, const Eigen::MatrixXd distanceMatrix, const unsigned int changedCoord) {
    unsigned int particle = unsigned(changedCoord/m_numberOfDimensions);

    m_positions                = positions;
    m_distanceMatrix           = distanceMatrix;
    calculateF                  (particle);
    calculateH                  (particle);
    calculateG                  (particle, changedCoord);
    calculateProbabilityRatio   (particle);
}

void PadeJastrow2::setArrays() {
    m_positionsOld          = m_positions;
    m_distanceMatrixOld     = m_distanceMatrix;
    m_hOldOld               = m_hOld;
    m_hOld                  = m_h;
    m_fOld                  = m_f;
    m_fSqrdOld              = m_fSqrd;
    m_gOld                  = m_g;
    m_gSqrdOld              = m_gSqrd;
    m_probabilityRatioOld   = m_probabilityRatio;
}

void PadeJastrow2::resetArrays() {
    m_positions             = m_positionsOld;
    m_distanceMatrix        = m_distanceMatrixOld;
    m_fSqrd                 = m_fSqrdOld;
    m_f                     = m_fOld;
    m_g                     = m_gOld;
    m_gSqrd                 = m_gSqrdOld;
    m_h                     = m_hOld;
    m_hOld                  = m_hOldOld;
    m_probabilityRatio      = m_probabilityRatioOld;
}

void PadeJastrow2::updateParameters(const Eigen::MatrixXd parameters, const unsigned short elementNumber) {
    m_elementNumber                     = elementNumber;
    m_maxNumberOfParametersPerElement   = m_system->getMaxNumberOfParametersPerElement();
    m_gamma                             = parameters(m_elementNumber, 0);
}

double PadeJastrow2::evaluateRatio() {
    return m_probabilityRatio;
}

double PadeJastrow2::computeGradient(const unsigned int k) {
    unsigned int k_p = unsigned(k/m_numberOfDimensions);  //Particle associated with k
    unsigned int k_d = k%m_numberOfDimensions;            //Dimension associated with k

    double derivative = 0;
    for(unsigned int i_p=0; i_p<k_p; i_p++) {
        unsigned int i = i_p * m_numberOfDimensions + k_d;
        derivative -= m_beta(i_p,k_p) * m_fSqrd(i_p,k_p) * m_g(i,k);
    }
    for(unsigned int j_p=k_p+1; j_p<m_numberOfParticles; j_p++) {
        unsigned int j = j_p * m_numberOfDimensions + k_d;
        derivative += m_beta(k_p,j_p) * m_fSqrd(k_p,j_p) * m_g(k,j);
    }
    return derivative;
}

double PadeJastrow2::computeLaplacian() {
    double derivative = 0;
    for(unsigned int k=0; k<m_numberOfFreeDimensions; k++) {
        unsigned int k_p = unsigned(k/m_numberOfDimensions);  //Particle associated with k
        unsigned int k_d = k%m_numberOfDimensions;            //Dimension associated with k
        for(unsigned int i_p=0; i_p<k_p; i_p++) {
            unsigned int i = i_p * m_numberOfDimensions + k_d;
            derivative -= m_beta(i_p,k_p) * m_fSqrd(i_p,k_p) * (1-(1+2*m_gamma*m_h(i_p,k_p))*m_gSqrd(i,k)) / m_distanceMatrix(i_p,k_p);
        }
        for(unsigned int j_p=k_p+1; j_p<m_numberOfParticles; j_p++) {
            unsigned int j = j_p * m_numberOfDimensions + k_d;
            derivative += m_beta(k_p,j_p) * m_fSqrd(k_p,j_p) * (1-(1+2*m_gamma*m_h(k_p,j_p))*m_gSqrd(k,j)) / m_distanceMatrix(k_p,j_p);
        }
    }
    return derivative;
}

Eigen::VectorXd PadeJastrow2::computeParameterGradient() {
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);

    double derivative = 0;
    for(unsigned int i_p=0; i_p<m_numberOfParticles; i_p++) {
        for(unsigned int j_p=i_p+1; j_p<m_numberOfParticles; j_p++) {
            derivative -= m_beta(i_p,j_p) * m_h(i_p, j_p) * m_h(i_p, j_p);
        }
    }
    gradients(0) = derivative;
    return gradients;
}
