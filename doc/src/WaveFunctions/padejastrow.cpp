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
    m_maxNumberOfParametersPerElement   = m_system->getMaxNumberOfParametersPerElement();
}


double PadeJastrow::calculateDistanceMatrixElement(const int i, const int j) {
    // Update element (i,j) in distance matrix

    double dist = 0;
    int parti   = m_numberOfDimensions*i;
    int partj   = m_numberOfDimensions*j;
    for(int d=0; d<m_numberOfDimensions; d++) {
        double diff = m_positions(parti+d)-m_positions(partj+d);
        dist += diff*diff;
    }
    return sqrt(dist);
}

void PadeJastrow::calculateDistanceMatrix() {
    m_distanceMatrix = Eigen::MatrixXd::Zero(m_numberOfParticles, m_numberOfParticles);
    for(int i=0; i<m_numberOfParticles; i++) {
        for(int j=i+1; j<m_numberOfParticles; j++) {
            m_distanceMatrix(i,j) = calculateDistanceMatrixElement(i,j);
            m_distanceMatrix(j,i) = m_distanceMatrix(i,j);
        }
    }
}

void PadeJastrow::calculateDistanceMatrixCross(const int par) {
    // Update distance matrix when position of particle "par" is changed
    for(int i=0; i<m_numberOfParticles; i++) {
        if(i!=par) {
            m_distanceMatrix(par, i) = calculateDistanceMatrixElement(par, i);
            m_distanceMatrix(i, par) = m_distanceMatrix(par, i);
        }
    }
}

void PadeJastrow::initializeBeta() {
    double a_sym, a_asym;

    if (m_numberOfDimensions == 2) {
        a_sym = 1. / 3;
        a_asym = 1.0;
    } else if (m_numberOfDimensions == 3) {
        a_sym = 1. / 4;
        a_asym = 1. / 2;
    } else {
        std::cout << "Unable to initialize Jastrow parameters: Unknown dimension" << std::endl;
        exit(0);
    }

    m_beta = Eigen::MatrixXd::Zero(m_numberOfParticles, m_numberOfParticles);

    int n2 = m_numberOfParticles/2;
    for (int i = 0; i < m_numberOfParticles; i++) {
        for (int j = 0; j < m_numberOfParticles; j++) {

            if ((j < n2 && i < n2) || (j >= n2 && i >= n2)) {
                m_beta(i, j) = a_sym;
            } else {
                m_beta(i, j) = a_asym;
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
}

void PadeJastrow::calculateG(int particle, int pRand) {
    //for(int j=0; j<m_numberOfFreeDimensions; j++) {
    //    m_g(pRand,j) = (m_positions(pRand) - m_positions(j))/m_distanceMatrix(particle,int(j/m_numberOfDimensions));
    //    m_g(j,pRand) = -m_g(pRand,j);
    //}
    for(int i=0; i<m_numberOfFreeDimensions; i++) {
        for(int j=0; j<m_numberOfFreeDimensions; j++) {
            m_g(i,j) = (m_positions(i) - m_positions(j))/m_distanceMatrix(int(i/m_numberOfDimensions),int(j/m_numberOfDimensions));
        }
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

void PadeJastrow::initializeArrays(const Eigen::VectorXd positions) {
    m_positions = positions;
    calculateDistanceMatrix();
    m_f     = (Eigen::MatrixXd::Ones(m_numberOfParticles, m_numberOfParticles) + m_gamma * m_distanceMatrix).cwiseInverse();
    m_g     = Eigen::MatrixXd::Zero(m_numberOfFreeDimensions, m_numberOfFreeDimensions);
    for(int i=0; i<m_numberOfFreeDimensions; i++) {
        for(int j=0; j<m_numberOfFreeDimensions; j++) {
            m_g(i,j) = (m_positions(i) - m_positions(j))/m_distanceMatrix(int(i/m_numberOfDimensions),int(j/m_numberOfDimensions));
        }
    }
    m_h = m_distanceMatrix.cwiseProduct(m_f);
    m_fOld = m_f;
    m_gOld = m_g;
    m_hOld = m_h;
    m_hOldOld = m_h;

    m_probabilityRatio = 1;

    initializeBeta();
}

void PadeJastrow::updateArrays(const Eigen::VectorXd positions, const int pRand) {
    int particle = int(pRand/m_numberOfDimensions);

    m_positionsOld          = m_positions;
    m_distanceMatrixOld     = m_distanceMatrix;
    m_hOldOld               = m_hOld;
    m_hOld                  = m_h;
    m_fOld                  = m_f;
    m_gOld                  = m_g;
    m_probabilityRatioOld   = m_probabilityRatio;

    m_positions             = positions;
    calculateDistanceMatrixCross(particle);
    calculateF                  (particle);
    calculateH                  (particle);
    calculateG                  (particle, pRand);
    calculateProbabilityRatio   (particle);
}

void PadeJastrow::resetArrays(int pRand) {
    m_positions         = m_positionsOld;
    m_distanceMatrix    = m_distanceMatrixOld;
    m_probabilityRatio  = m_probabilityRatioOld;
    m_f                 = m_fOld;
    m_g                 = m_gOld;
    m_h                 = m_hOld;
    m_hOld              = m_hOldOld;
}

void PadeJastrow::updateParameters(const Eigen::MatrixXd parameters, const int elementNumber) {
    m_elementNumber = elementNumber;
    m_gamma         = parameters(m_elementNumber, 0);
}

double PadeJastrow::evaluateRatio() {
    return m_probabilityRatio;
}

double PadeJastrow::computeFirstDerivative(const int k) {
    int k_p = int(k/m_numberOfDimensions);  //Particle associated with k
    int k_d = k%m_numberOfDimensions;       //Dimension associated with k

    double derivative = 0;
    for(int j_p=0; j_p<m_numberOfParticles; j_p++) {
        int j = j_p * m_numberOfDimensions + k_d;
        if(j_p!=k_p) {
            derivative += m_beta(k_p,j_p) * m_f(k_p,j_p) * m_f(k_p, j_p) * m_g(k,j);
        }
    }
    return derivative;
}

double PadeJastrow::computeSecondDerivative() {
    double derivative = 0;
    for(int i=0; i<m_numberOfFreeDimensions; i++) {
        int i_p = int(i/m_numberOfDimensions);  //Particle associated with k
        int i_d = i%m_numberOfDimensions;       //Dimension associated with k
        for(int j_p=0; j_p<m_numberOfParticles; j_p++) {
            int j = j_p * m_numberOfDimensions + i_d;
            if(j_p!=i_p) {
                derivative += m_beta(i_p,j_p) * m_f(i_p,j_p) * m_f(i_p,j_p) * (1-(1+2*m_gamma*m_h(i_p,j_p))*m_g(i,j)*m_g(i,j)) / m_distanceMatrix(i_p,j_p);
            }
        }
    }
    return derivative;
}

Eigen::VectorXd PadeJastrow::computeFirstEnergyDerivative(const int k) {
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);

    int k_p = int(k/m_numberOfDimensions);  //Particle associated with k
    int k_d = k%m_numberOfDimensions;       //Dimension associated with k

    double derivative = 0;
    for(int j_p=0; j_p<m_numberOfParticles; j_p++) {
        int j = j_p * m_numberOfDimensions + k_d;
        derivative += m_beta(k_p,j_p) * m_f(k_p, j_p) * m_f(k_p, j_p) * m_f(k_p, j_p) * (m_positions(k) - m_positions(j));
    }
    gradients(0) = derivative;
    return gradients;
}

Eigen::VectorXd PadeJastrow::computeSecondEnergyDerivative() {
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);

    double derivative = 0;
    for(int i=0; i<m_numberOfFreeDimensions; i++) {
        int i_p = int(i/m_numberOfDimensions);  //Particle associated with k
        int i_d = i%m_numberOfDimensions;       //Dimension associated with k
        for(int j_p=0; j_p<m_numberOfParticles; j_p++) {
            int j = j_p * m_numberOfDimensions + i_d;
            if(j_p!=i_p) {
                derivative += m_beta(i_p,j_p) * m_f(i_p,j_p) * m_f(i_p,j_p) * m_f(i_p,j_p) * (1 - 4 * m_gamma * m_h(i_p,j_p) * m_g(i,j) * m_g(i,j));
            }
        }
    }
    gradients(0) = derivative;
    return gradients;
}
