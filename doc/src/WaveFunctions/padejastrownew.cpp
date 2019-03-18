#include "padejastrownew.h"
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include <iostream>

PadeJastrowNew::PadeJastrowNew(System* system) :
        WaveFunction(system) {
    m_numberOfParticles                 = m_system->getNumberOfParticles();
    m_numberOfDimensions                = m_system->getNumberOfDimensions();
    m_numberOfFreeDimensions            = m_system->getNumberOfFreeDimensions();
    m_maxNumberOfParametersPerElement   = m_system->getMaxNumberOfParametersPerElement();
}

double PadeJastrowNew::calculateDistanceMatrixElement(const int i, const int j, const Eigen::VectorXd particles) {
    // Update element (i,j) in Dist_inv_invance matrix

    double dist = 0;
    for(int d=0; d<m_numberOfDimensions; d++) {
        double diff = particles(m_numberOfDimensions*i+d)-particles(m_numberOfDimensions*j+d);
        dist += diff*diff;
    }
    return sqrt(dist);
}

Eigen::MatrixXd PadeJastrowNew::calculateDistanceMatrix(const Eigen::VectorXd particles) {
    Eigen::MatrixXd distanceMatrix = Eigen::MatrixXd::Zero(m_numberOfParticles, m_numberOfParticles);
    for(int i=0; i<m_numberOfParticles; i++) {
        for(int j=0; j<m_numberOfParticles; j++) {
            if(i!=j) {
                distanceMatrix(i,j) = calculateDistanceMatrixElement(i,j,particles);
            }
        }
    }
    return distanceMatrix;
}

void PadeJastrowNew::calculateDistanceMatrixCross(const int par, const Eigen::VectorXd particles, Eigen::MatrixXd &distanceMatrix) {
    // Update distance matrix when position of particle "par" is changed
    for(int i=0; i<m_numberOfParticles; i++) {
        if(i!=par) {
            distanceMatrix(par, i) = calculateDistanceMatrixElement(par, i, particles);
            distanceMatrix(i, par) = calculateDistanceMatrixElement(i, par, particles);
        }
    }
}

void PadeJastrowNew::initializeArrays(const Eigen::VectorXd positions) {
    m_positions = positions;
    m_distanceMatrix = calculateDistanceMatrix(positions);
    m_f     = (Eigen::MatrixXd::Ones(m_numberOfParticles, m_numberOfParticles) + m_gamma * m_distanceMatrix).cwiseInverse();
    m_g     = Eigen::MatrixXd::Zero(m_numberOfFreeDimensions, m_numberOfFreeDimensions);
    for(int i=0; i<m_numberOfFreeDimensions; i++) {
        for(int j=0; j<m_numberOfFreeDimensions; j++) {
            m_g(i,j) = (m_positions(i) - m_positions(j))/m_distanceMatrix(int(i/m_numberOfDimensions),int(j/m_numberOfDimensions));
        }
    }
    m_h = m_distanceMatrix.cwiseProduct(m_f);
    m_hOld = m_h;

    m_ratio = 1;

    double a_sym, a_asym;

    if (m_numberOfDimensions == 2) {
        a_sym = 1. / 3;
        a_asym = 1.0;
    } else if (m_numberOfDimensions == 3) {
        a_sym = 1. / 4;
        a_asym = 1. / 2;
    } else {
        std::cout << "Unable to initialize Jastrow paremters: Unknown dimension" << std::endl;
        exit(1);
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

void PadeJastrowNew::updateArrays(const Eigen::VectorXd positions, const int pRand) {
    int particle = int(pRand/m_numberOfDimensions);

    m_oldPositions = m_positions;
    m_positions = positions;

    m_oldDistanceMatrix = m_distanceMatrix;
    calculateDistanceMatrixCross(particle, positions, m_distanceMatrix);

    m_f     = (Eigen::MatrixXd::Ones(m_numberOfParticles, m_numberOfParticles) + m_gamma * m_distanceMatrix).cwiseInverse();

    for(int i=0; i<m_numberOfFreeDimensions; i++) {
        for(int j=0; j<m_numberOfFreeDimensions; j++) {
            m_g(i,j) = (m_positions(i) - m_positions(j))/m_distanceMatrix(int(i/m_numberOfDimensions),int(j/m_numberOfDimensions));
        }
    }

    m_hOld = m_h;
    m_h = m_distanceMatrix.cwiseProduct(m_f);

    m_oldRatio = m_ratio;
    double ratio = 0;
    for(int i=particle; i<m_numberOfParticles; i++) {
        ratio += m_beta(particle, i) * (m_h(particle, i) - m_hOld(particle, i));
    }
    m_ratio = exp(2*ratio);
}

void PadeJastrowNew::resetArrays(int pRand) {
    m_positions      = m_oldPositions;
    m_distanceMatrix = m_oldDistanceMatrix;
    m_ratio          = m_oldRatio;
}

void PadeJastrowNew::updateParameters(const Eigen::MatrixXd parameters, const int elementNumber) {
    m_elementNumber = elementNumber;
    m_gamma = parameters(m_elementNumber, 0);
}

double PadeJastrowNew::evaluateRatio() {
    return m_ratio;
}

double PadeJastrowNew::computeFirstDerivative(const int k) {
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

double PadeJastrowNew::computeSecondDerivative() {
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

Eigen::VectorXd PadeJastrowNew::computeFirstEnergyDerivative(const int k) {
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

Eigen::VectorXd PadeJastrowNew::computeSecondEnergyDerivative() {
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
