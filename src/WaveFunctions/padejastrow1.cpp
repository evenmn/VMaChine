#include "padejastrow.h"
#include "../system.h"
#include "wavefunction.h"
#include <cassert>
#include <iostream>

PadeJastrow::PadeJastrow(System *system)
    : WaveFunction(system)
{
    m_numberOfParticles = m_system->getNumberOfParticles();
    m_numberOfDimensions = m_system->getNumberOfDimensions();
    m_degreesOfFreedom = m_system->getNumberOfFreeDimensions();
}

void PadeJastrow::setConstants(const int elementNumber)
{
    m_elementNumber = elementNumber;
    m_gradients = Eigen::VectorXd::Zero(m_system->getMaxParameters());
}

void PadeJastrow::initializeArrays(const Eigen::VectorXd positions,
                                   const Eigen::VectorXd radialVector,
                                   const Eigen::MatrixXd distanceMatrix)
{
    m_positions = positions;
    m_distanceMatrix = distanceMatrix;
    m_probabilityRatio = 1;
    m_f = (Eigen::MatrixXd::Ones(m_numberOfParticles, m_numberOfParticles)
           + m_gamma * m_distanceMatrix)
              .cwiseInverse();
    m_h = m_distanceMatrix.cwiseProduct(m_f);
    m_fSqrd = m_f.cwiseAbs2();
    initializePrincipalDistance();
    initializeBeta();
}

void PadeJastrow::updateArrays(const Eigen::VectorXd positions,
                               const Eigen::VectorXd radialVector,
                               const Eigen::MatrixXd distanceMatrix,
                               const int changedCoord)
{
    int particle = int(changedCoord / m_numberOfDimensions);
    m_positions = positions;
    m_distanceMatrix = distanceMatrix;
    calculateF(particle);
    calculateH(particle);
    updatePrincipalDistance(changedCoord, particle);
    calculateProbabilityRatio(particle);
}

void PadeJastrow::setArrays()
{
    m_positionsOld = m_positions;
    m_distanceMatrixOld = m_distanceMatrix;
    m_hOldOld = m_hOld;
    m_hOld = m_h;
    m_fOld = m_f;
    m_fSqrdOld = m_fSqrd;
    m_principalDistanceOld = m_principalDistance;
    m_probabilityRatioOld = m_probabilityRatio;
}

void PadeJastrow::resetArrays()
{
    m_positions = m_positionsOld;
    m_distanceMatrix = m_distanceMatrixOld;
    m_fSqrd = m_fSqrdOld;
    m_f = m_fOld;
    m_principalDistance = m_principalDistanceOld;
    m_h = m_hOld;
    m_hOld = m_hOldOld;
    m_probabilityRatio = m_probabilityRatioOld;
}

void PadeJastrow::updateParameters(const Eigen::MatrixXd parameters)
{
    m_gamma = parameters(m_elementNumber, 0);
}

double PadeJastrow::evaluateRatio()
{
    return m_probabilityRatio;
}

double PadeJastrow::computeGradient(const int k)
{
    int k_p = int(k / m_numberOfDimensions); //Particle associated with k
    int k_d = k % m_numberOfDimensions;      //Dimension associated with k

    double derivative = 0;
    for (int j_p = 0; j_p < m_numberOfParticles; j_p++) {
        int j = j_p * m_numberOfDimensions + k_d;
        if (j_p != k_p) {
            derivative += m_beta(k_p, j_p) * m_fSqrd(k_p, j_p) * m_principalDistance(k, j);
        }
    }
    return derivative;
}

double PadeJastrow::computeLaplacian()
{
    double derivative = 0;
    for (int i = 0; i < m_degreesOfFreedom; i++) {
        int i_p = int(i / m_numberOfDimensions); //Particle associated with k
        int i_d = i % m_numberOfDimensions;      //Dimension associated with k
        for (int j_p = i_p + 1; j_p < m_numberOfParticles; j_p++) {
            int j = j_p * m_numberOfDimensions + i_d;
            derivative += m_beta(i_p, j_p) * m_fSqrd(i_p, j_p)
                          * (1
                             - (1 + 2 * m_gamma * m_h(i_p, j_p)) * m_principalDistance(i, j)
                                   * m_principalDistance(i, j))
                          / m_distanceMatrix(i_p, j_p);
        }
    }
    return 2 * derivative;
}

Eigen::VectorXd PadeJastrow::computeParameterGradient()
{
    //double derivative2 = m_beta.cwiseProduct(m_h.cwiseAbs2()).sum();
    double derivative = 0;
    for (int i_p = 0; i_p < m_numberOfParticles; i_p++) {
        for (int j_p = i_p + 1; j_p < m_numberOfParticles; j_p++) {
            derivative -= m_beta(i_p, j_p) * m_h(i_p, j_p) * m_h(i_p, j_p);
        }
    }
    m_gradients(0) = derivative;
    return m_gradients;
}

void PadeJastrow::initializePrincipalDistance()
{
    m_principalDistance = Eigen::MatrixXd::Zero(m_degreesOfFreedom, m_degreesOfFreedom);
    for (int n = 1; n < m_numberOfParticles; n++) {
        int step = n * m_numberOfDimensions;
        for (int j = 0; j < m_degreesOfFreedom - step; j++) {
            int p = int(j / m_numberOfDimensions);
            m_principalDistance(j, j + step) = (m_positions(j) - m_positions(j + step))
                                               / m_distanceMatrix(p, p + n);
            m_principalDistance(j + step, j) = -m_principalDistance(j, j + step);
        }
    }
}

void PadeJastrow::updatePrincipalDistance(int i, int i_p)
{
    /* Update of the principal distance matrix
     * Arguments:
     * 
     * {int} i:     The changed coordinate
     * {int} i_p:   The moved particle
     */
    int i_d = i % m_numberOfDimensions;
    for (int j_p = 0; j_p < i_p; j_p++) {
        int j = i_d + j_p * m_numberOfDimensions;
        m_principalDistance(i, j) = (m_positions(i) - m_positions(j)) / m_distanceMatrix(i_p, j_p);
        m_principalDistance(j, i) = -m_principalDistance(i, j);
    }
    for (int j_p = i_p + 1; j_p < m_numberOfParticles; j_p++) {
        int j = i_d + j_p * m_numberOfDimensions;
        m_principalDistance(i, j) = (m_positions(i) - m_positions(j)) / m_distanceMatrix(i_p, j_p);
        m_principalDistance(j, i) = -m_principalDistance(i, j);
    }
}

void PadeJastrow::initializeBeta()
{
    double symmetric, antisymmetric;
    if (m_numberOfDimensions == 2) {
        symmetric = 1. / 3;
        antisymmetric = 1. / 1;
    } else if (m_numberOfDimensions == 3) {
        symmetric = 1. / 4;
        antisymmetric = 1. / 2;
    } else if (m_numberOfDimensions == 4) {
        symmetric = 1. / 5;
        antisymmetric = 1. / 3;
    } else {
        std::cout << "Number of dimensions should be in the range [1, 4]" << std::endl;
        exit(0);
    }

    m_beta = Eigen::MatrixXd::Zero(m_numberOfParticles, m_numberOfParticles);

    int numberOfParticlesHalf = m_numberOfParticles / 2;
    for (int i = 0; i < m_numberOfParticles; i++) {
        for (int j = 0; j < m_numberOfParticles; j++) {
            if ((j < numberOfParticlesHalf && i < numberOfParticlesHalf)
                || (j >= numberOfParticlesHalf && i >= numberOfParticlesHalf)) {
                m_beta(i, j) = symmetric;
            } else {
                m_beta(i, j) = antisymmetric;
            }
        }
    }
}

void PadeJastrow::calculateF(int i_p)
{
    //m_f     = (Eigen::MatrixXd::Ones(m_numberOfParticles, m_numberOfParticles) + m_gamma * m_distanceMatrix).cwiseInverse();

    for (int i = 0; i < m_numberOfParticles; i++) {
        m_f(i, i_p) = 1 / (1 + m_gamma * m_distanceMatrix(i, i_p));
        m_f(i_p, i) = m_f(i, i_p);
    }
    m_fSqrd = m_f.cwiseAbs2();
}

void PadeJastrow::calculateH(int i_p)
{
    //m_h = m_distanceMatrix.cwiseProduct(m_f);

    for (int i = 0; i < m_numberOfParticles; i++) {
        m_h(i_p, i) = m_distanceMatrix(i_p, i) * m_f(i_p, i);
        m_h(i, i_p) = m_h(i_p, i);
    }
}

void PadeJastrow::calculateProbabilityRatio(int i_p)
{
    //double ratio = double(m_beta.row(particle).transpose() * (m_h.row(particle) - m_hOld.row(particle)));

    double ratio = 0;
    for (int i = i_p; i < m_numberOfParticles; i++) {
        ratio += m_beta(i_p, i) * (m_h(i_p, i) - m_hOld(i_p, i));
    }
    m_probabilityRatio = exp(2 * ratio);
}
