#include "padejastrow.h"
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include <iostream>

PadeJastrow::PadeJastrow(System* system) :
        WaveFunction(system) {
    m_numberOfParticles                 = m_system->getNumberOfParticles();
    m_numberOfDimensions                = m_system->getNumberOfDimensions();
    m_numberOfFreeDimensions            = m_system->getNumberOfFreeDimensions();
}

void PadeJastrow::initializeBeta() {
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

    int numberOfParticlesHalf = m_numberOfParticles/2;
    for (int i = 0; i < m_numberOfParticles; i++) {
        for (int j = 0; j < m_numberOfParticles; j++) {

            if ((j <  numberOfParticlesHalf && i <  numberOfParticlesHalf) ||  \
                (j >= numberOfParticlesHalf && i >= numberOfParticlesHalf)) {
                m_beta(i, j) = symmetric;
            } else {
                m_beta(i, j) = antisymmetric;
            }
        }
    }
}

void PadeJastrow::calculateF(int particle) {
    //m_f     = (Eigen::MatrixXd::Ones(m_numberOfParticles, m_numberOfParticles) + m_gamma * m_distanceMatrix).cwiseInverse();

    for(int i=0; i<m_numberOfParticles; i++) {
        m_f(i, particle) = 1/(1 + m_gamma * m_distanceMatrix(i, particle));
        m_f(particle, i) = m_f(i, particle);
    }
    m_fSqrd = m_f.cwiseAbs2();
}

void PadeJastrow::calculateG(int pRand) {
    for(int i=0; i<m_numberOfFreeDimensions; i++) {
        m_g(pRand,i) = m_positions(pRand) - m_positions(i);
        m_g(i,pRand) = -m_g(pRand,i);
    }
}

void PadeJastrow::calculateH(int particle) {
    //m_h = m_distanceMatrix.cwiseProduct(m_f);

    for(int i=0; i<m_numberOfParticles; i++) {
        m_h(particle, i) = m_distanceMatrix(particle, i) * m_f(particle, i);
        m_h(i, particle) = m_h(particle, i);
    }
}

void PadeJastrow::calculateProbabilityRatio(int particle) {
    //double ratio = (m_beta.row(particle).transpose() * (m_h.row(particle) - m_hOld.row(particle))).sum();

    double ratio = 0;
    for(int i=particle; i<m_numberOfParticles; i++) {
        ratio += m_beta(particle, i) * (m_h(particle, i) - m_hOld(particle, i));
    }
    m_probabilityRatio = exp(2*ratio);
}

void PadeJastrow::initializeArrays(const Eigen::VectorXd positions, const Eigen::VectorXd radialVector, const Eigen::MatrixXd distanceMatrix) {
    m_positions         = positions;
    m_distanceMatrix    = distanceMatrix;
    m_probabilityRatio  = 1;
    m_f     = (Eigen::MatrixXd::Ones(m_numberOfParticles, m_numberOfParticles) + m_gamma * m_distanceMatrix).cwiseInverse();
    m_g     = Eigen::MatrixXd::Zero(m_numberOfFreeDimensions, m_numberOfFreeDimensions);
    for(int i=0; i<m_numberOfFreeDimensions; i++) {
        for(int j=i; j<m_numberOfFreeDimensions; j++) {
            m_g(i,j) = m_positions(i) - m_positions(j);
            m_g(j,i) = -m_g(i,j);
        }
    }
    m_h = m_distanceMatrix.cwiseProduct(m_f);
    m_fSqrd = m_f.cwiseAbs2();

    setArrays();
    initializeBeta();
}

void PadeJastrow::updateArrays(const Eigen::VectorXd positions, const Eigen::VectorXd radialVector, const Eigen::MatrixXd distanceMatrix, const int changedCoord) {
    int particle = int(changedCoord/m_numberOfDimensions);
    setArrays();

    m_positions                = positions;
    m_distanceMatrix           = distanceMatrix;
    calculateF                  (particle);
    calculateH                  (particle);
    calculateG                  (changedCoord);
    calculateProbabilityRatio   (particle);
}

void PadeJastrow::setArrays() {
    m_positionsOld          = m_positions;
    m_distanceMatrixOld     = m_distanceMatrix;
    m_hOldOld               = m_hOld;
    m_hOld                  = m_h;
    m_fOld                  = m_f;
    m_fSqrdOld              = m_fSqrd;
    m_gOld                  = m_g;
    m_probabilityRatioOld   = m_probabilityRatio;
}

void PadeJastrow::resetArrays() {
    m_positions             = m_positionsOld;
    m_distanceMatrix        = m_distanceMatrixOld;
    m_fSqrd                 = m_fSqrdOld;
    m_f                     = m_fOld;
    m_g                     = m_gOld;
    m_h                     = m_hOld;
    m_hOld                  = m_hOldOld;
    m_probabilityRatio      = m_probabilityRatioOld;
}

void PadeJastrow::updateParameters(const Eigen::MatrixXd parameters, const int elementNumber) {
    m_elementNumber                     = elementNumber;
    m_maxNumberOfParametersPerElement   = m_system->getMaxNumberOfParametersPerElement();
    m_gamma                             = parameters(m_elementNumber, 0);
}

double PadeJastrow::evaluateRatio() {
    return m_probabilityRatio;
}

double PadeJastrow::computeGradient(const int k) {
    int k_p = int(k/m_numberOfDimensions);  //Particle associated with k
    int k_d = k%m_numberOfDimensions;       //Dimension associated with k

    double derivative = 0;
    for(int j_p=0; j_p<m_numberOfParticles; j_p++) {
        int j = j_p * m_numberOfDimensions + k_d;
        if(j_p!=k_p) {
            derivative += m_beta(k_p,j_p) * m_fSqrd(k_p,j_p) * m_g(k,j)/m_distanceMatrix(k_p,j_p);
        }
    }
    return derivative;
}

double PadeJastrow::computeLaplacian() {
    double derivative = 0;
    for(int i=0; i<m_numberOfFreeDimensions; i++) {
        int i_p = int(i/m_numberOfDimensions);  //Particle associated with k
        int i_d = i%m_numberOfDimensions;       //Dimension associated with k
        for(int j_p=i_p+1; j_p<m_numberOfParticles; j_p++) {
            int j = j_p * m_numberOfDimensions + i_d;
            derivative += m_beta(i_p,j_p) * m_fSqrd(i_p,j_p) * (1-(1+2*m_gamma*m_h(i_p,j_p))*m_g(i,j)*m_g(i,j) / (m_distanceMatrix(i_p,j_p)*m_distanceMatrix(i_p,j_p))) / m_distanceMatrix(i_p,j_p);
        }
    }
    return 2 * derivative;
}

Eigen::VectorXd PadeJastrow::computeParameterGradient() {
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);

    double derivative = 0;
    for(int i_p=0; i_p<m_numberOfParticles; i_p++) {
        for(int j_p=i_p+1; j_p<m_numberOfParticles; j_p++) {
            derivative -= m_beta(i_p,j_p) * m_h(i_p, j_p) * m_h(i_p, j_p);
        }
    }
    gradients(0) = derivative;
    return gradients;
}
