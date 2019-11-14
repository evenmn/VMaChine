#include "padejastrow.h"
#include "../system.h"
#include "wavefunction.h"
#include <cassert>
#include <iostream>

PadeJastrow::PadeJastrow(System *system)
    : WaveFunction(system)
{}

void PadeJastrow::setConstants(const int elementNumber)
{
    m_elementNumber = elementNumber;
    m_numberOfParticles = m_system->getNumberOfParticles();
    m_numberOfDimensions = m_system->getNumberOfDimensions();
    m_degreesOfFreedom = m_system->getNumberOfFreeDimensions();
}

void PadeJastrow::initializeArrays(const Eigen::VectorXd positions,
                                   const Eigen::VectorXd /*radialVector*/,
                                   const Eigen::MatrixXd distanceMatrix)
{
    m_positions = positions;
    m_distanceMatrix = distanceMatrix;
    m_probabilityRatio = 1;
    initializePrincipalDistance();
    initializeBeta();
    initializeMatrices();
}

void PadeJastrow::updateArrays(const Eigen::VectorXd positions,
                               const Eigen::VectorXd /*radialVector*/,
                               const Eigen::MatrixXd distanceMatrix,
                               const int changedCoord)
{
    int particle = int(changedCoord / m_numberOfDimensions);

    m_positions = positions;
    m_distanceMatrix = distanceMatrix;
    updateMatrices(particle);
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
    m_principalDistanceOld = m_principalDistance;
    m_probabilityRatioOld = m_probabilityRatio;
}

void PadeJastrow::resetArrays()
{
    m_positions = m_positionsOld;
    m_distanceMatrix = m_distanceMatrixOld;
    m_f = m_fOld;
    m_h = m_hOld;
    m_hOld = m_hOldOld;
    m_principalDistance = m_principalDistanceOld;
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
    int k_p = int(k / m_numberOfDimensions);
    int k_d = k % m_numberOfDimensions;

    double derivative = 0;
    for (int j_p = 0; j_p < k_p; j_p++) {
        int j = j_p * m_numberOfDimensions + k_d;
        derivative += m_f(k_p, j_p) * m_principalDistance(k, j);
    }
    for (int j_p = k_p + 1; j_p < m_numberOfParticles; j_p++) {
        int j = j_p * m_numberOfDimensions + k_d;
        derivative += m_f(k_p, j_p) * m_principalDistance(k, j);
    }
    return derivative;
}

double PadeJastrow::computeLaplacian()
{
    double derivative = 0;
    for (int k = 0; k < m_degreesOfFreedom; k++) {
        int k_p = int(k / m_numberOfDimensions);
        int k_d = k % m_numberOfDimensions;
        for (int j_p = k_p + 1; j_p < m_numberOfParticles; j_p++) {
            int j = j_p * m_numberOfDimensions + k_d;

            derivative += m_f(k_p, j_p)
                          * (1
                             - (1 + 2 * m_gamma * m_h(k_p, j_p)) * m_principalDistance(k, j)
                                   * m_principalDistance(k, j))
                          / m_distanceMatrix(k_p, j_p);
        }
    }
    return 2 * derivative;
}

Eigen::VectorXd PadeJastrow::computeParameterGradient()
{
    m_gradients = Eigen::VectorXd::Zero(m_system->getMaxParameters());
    m_gradients(0) = -m_beta.cwiseProduct(m_h.cwiseAbs2()).sum() / 2;
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

void PadeJastrow::initializeMatrices()
{
    Eigen::MatrixXd f = (Eigen::MatrixXd::Ones(m_numberOfParticles, m_numberOfParticles)
                         + m_gamma * m_distanceMatrix)
                            .cwiseInverse();
    m_h = m_distanceMatrix.cwiseProduct(f);
    m_f = m_beta.cwiseProduct(f.cwiseAbs2());
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

void PadeJastrow::updatePrincipalDistance(int i, int i_p)
{
    /* Update of the principal distance matrix
     * Arguments:
     * 
     * {int} i:     The changed coordinate
     * {int} i_p:   The moved particle
     */

    /*
    for (int d = 0; d < m_numberOfDimensions; d++) {
        int i = i_p + d;
        for (int j_p = 1; j_p < i_p; j_p++) {
            int j = d + j_p * m_numberOfDimensions;
            m_principalDistance(i, j) = (m_positions(i) - m_positions(j))
                                        / m_distanceMatrix(i_p, j_p);
            m_principalDistance(j, i) = -m_principalDistance(i, j);
        }
        for (int j_p = i_p + 1; j_p < m_numberOfParticles; j_p++) {
            int j = d + j_p * m_numberOfDimensions;
            m_principalDistance(i, j) = (m_positions(i) - m_positions(j))
                                        / m_distanceMatrix(i_p, j_p);
            m_principalDistance(j, i) = -m_principalDistance(i, j);
        }
    }
    */

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

void PadeJastrow::updateMatrices(int i_p)
{
    for (int j_p = 0; j_p < m_numberOfParticles; j_p++) {
        double f = 1 / (1 + m_gamma * m_distanceMatrix(i_p, j_p));
        m_h(i_p, j_p) = m_distanceMatrix(i_p, j_p) * f;
        m_h(j_p, i_p) = m_h(i_p, j_p);
        m_f(i_p, j_p) = m_beta(i_p, j_p) * f * f;
        m_f(j_p, i_p) = m_f(i_p, j_p);
    }
}

void PadeJastrow::calculateProbabilityRatio(int i_p)
{
    double ratio = 0;
    for (int i = i_p; i < m_numberOfParticles; i++) {
        ratio += m_beta(i_p, i) * (m_h(i_p, i) - m_hOld(i_p, i));
    }
    m_probabilityRatio = exp(2 * ratio);
}
