#include "padejastrow.h"
#include "../system.h"
#include "wavefunction.h"
#include <cassert>
#include <iostream>

PadeJastrow::PadeJastrow(System *system)
    : WaveFunction(system)
{}

void PadeJastrow::setConstants(const arma::uword elementNumber)
{
    m_elementNumber = elementNumber;
    m_numberOfParticles = m_system->getNumberOfParticles();
    m_numberOfDimensions = m_system->getNumberOfDimensions();
    m_degreesOfFreedom = m_system->getNumberOfFreeDimensions();
    m_gradients.zeros(m_system->getMaxParameters());
}

void PadeJastrow::initializeArrays(const arma::vec positions,
                                   const arma::vec radialVector,
                                   const arma::mat distanceMatrix)
{
    m_positions = positions;
    m_distanceMatrix = distanceMatrix;
    m_probabilityRatio = 1;
    initializePrincipalDistance();
    initializeBeta();
    initializeMatrices();
}

void PadeJastrow::updateArrays(const arma::vec positions,
                               const arma::vec radialVector,
                               const arma::mat distanceMatrix,
                               const arma::uword changedCoord)
{
    arma::uword particle = arma::uword(changedCoord / m_numberOfDimensions);

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

void PadeJastrow::updateParameters(const arma::mat parameters)
{
    m_gamma = parameters(m_elementNumber, 0);
}

double PadeJastrow::evaluateRatio()
{
    return m_probabilityRatio;
}

double PadeJastrow::computeGradient(const arma::uword k)
{
    arma::uword k_p = arma::uword(k / m_numberOfDimensions);
    arma::uword k_d = k % m_numberOfDimensions;

    double derivative = 0;
    for (arma::uword j_p = 0; j_p < k_p; j_p++) {
        arma::uword j = j_p * m_numberOfDimensions + k_d;
        derivative += m_f(k_p, j_p) * m_principalDistance(k, j);
    }
    for (arma::uword j_p = k_p + 1; j_p < m_numberOfParticles; j_p++) {
        arma::uword j = j_p * m_numberOfDimensions + k_d;
        derivative += m_f(k_p, j_p) * m_principalDistance(k, j);
    }
    return derivative;
}

double PadeJastrow::computeLaplacian()
{
    double derivative = 0;
    for (arma::uword k = 0; k < m_degreesOfFreedom; k++) {
        arma::uword k_p = arma::uword(k / m_numberOfDimensions);
        arma::uword k_d = k % m_numberOfDimensions;
        for (arma::uword j_p = k_p + 1; j_p < m_numberOfParticles; j_p++) {
            arma::uword j = j_p * m_numberOfDimensions + k_d;

            derivative += m_f(k_p, j_p)
                          * (1
                             - (1 + 2 * m_gamma * m_h(k_p, j_p)) * m_principalDistance(k, j)
                                   * m_principalDistance(k, j))
                          / m_distanceMatrix(k_p, j_p);
        }
    }
    return 2 * derivative;
}

arma::vec PadeJastrow::computeParameterGradient()
{
    m_gradients(0) = - arma::accu(m_beta % arma::square(m_h)) / 2;
    return m_gradients;
}

void PadeJastrow::initializePrincipalDistance()
{
    m_principalDistance.zeros(m_degreesOfFreedom, m_degreesOfFreedom);
    for (arma::uword n = 1; n < m_numberOfParticles; n++) {
        arma::uword step = n * m_numberOfDimensions;
        for (arma::uword j = 0; j < m_degreesOfFreedom - step; j++) {
            arma::uword p = arma::uword(j / m_numberOfDimensions);
            m_principalDistance(j, j + step) = (m_positions(j) - m_positions(j + step))
                                               / m_distanceMatrix(p, p + n);
            m_principalDistance(j + step, j) = -m_principalDistance(j, j + step);
        }
    }
}

void PadeJastrow::initializeMatrices()
{
    arma::mat f = arma::pow(arma::ones<arma::mat>(m_numberOfParticles, m_numberOfParticles) + m_gamma * m_distanceMatrix, -1);
    m_h = m_distanceMatrix % f;
    m_f = m_beta % arma::square(f);
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

    m_beta.zeros(m_numberOfParticles, m_numberOfParticles);

    arma::uword numberOfParticlesHalf = m_numberOfParticles / 2;
    for (arma::uword i = 0; i < m_numberOfParticles; i++) {
        for (arma::uword j = 0; j < m_numberOfParticles; j++) {
            if ((j < numberOfParticlesHalf && i < numberOfParticlesHalf)
                || (j >= numberOfParticlesHalf && i >= numberOfParticlesHalf)) {
                m_beta(i, j) = symmetric;
            } else {
                m_beta(i, j) = antisymmetric;
            }
        }
    }
}

void PadeJastrow::updatePrincipalDistance(arma::uword i, arma::uword i_p)
{
    /* Update of the principal distance matrix
     * Arguments:
     * 
     * {arma::uword} i:     The changed coordinate
     * {arma::uword} i_p:   The moved particle
     */

    /*
    for (arma::uword d = 0; d < m_numberOfDimensions; d++) {
        arma::uword i = i_p + d;
        for (arma::uword j_p = 1; j_p < i_p; j_p++) {
            arma::uword j = d + j_p * m_numberOfDimensions;
            m_principalDistance(i, j) = (m_positions(i) - m_positions(j))
                                        / m_distanceMatrix(i_p, j_p);
            m_principalDistance(j, i) = -m_principalDistance(i, j);
        }
        for (arma::uword j_p = i_p + 1; j_p < m_numberOfParticles; j_p++) {
            arma::uword j = d + j_p * m_numberOfDimensions;
            m_principalDistance(i, j) = (m_positions(i) - m_positions(j))
                                        / m_distanceMatrix(i_p, j_p);
            m_principalDistance(j, i) = -m_principalDistance(i, j);
        }
    }
    */

    m_principalDistance.zeros(m_degreesOfFreedom, m_degreesOfFreedom);
    for (arma::uword n = 1; n < m_numberOfParticles; n++) {
        arma::uword step = n * m_numberOfDimensions;
        for (arma::uword j = 0; j < m_degreesOfFreedom - step; j++) {
            arma::uword p = arma::uword(j / m_numberOfDimensions);
            m_principalDistance(j, j + step) = (m_positions(j) - m_positions(j + step))
                                               / m_distanceMatrix(p, p + n);
            m_principalDistance(j + step, j) = -m_principalDistance(j, j + step);
        }
    }
}

void PadeJastrow::updateMatrices(arma::uword i_p)
{
    for (arma::uword j_p = 0; j_p < m_numberOfParticles; j_p++) {
        double f = 1 / (1 + m_gamma * m_distanceMatrix(i_p, j_p));
        m_h(i_p, j_p) = m_distanceMatrix(i_p, j_p) * f;
        m_h(j_p, i_p) = m_h(i_p, j_p);
        m_f(i_p, j_p) = m_beta(i_p, j_p) * f * f;
        m_f(j_p, i_p) = m_f(i_p, j_p);
    }
}

void PadeJastrow::calculateProbabilityRatio(arma::uword i_p)
{
    double ratio = 0;
    for (arma::uword i = i_p; i < m_numberOfParticles; i++) {
        ratio += m_beta(i_p, i) * (m_h(i_p, i) - m_hOld(i_p, i));
    }
    m_probabilityRatio = exp(2 * ratio);
}
