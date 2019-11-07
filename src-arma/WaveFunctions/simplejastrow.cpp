#include "simplejastrow.h"
#include "../system.h"
#include "wavefunction.h"
#include <cassert>
#include <iostream>

SimpleJastrow::SimpleJastrow(System *system)
    : WaveFunction(system)
{}

void SimpleJastrow::setConstants(const arma::uword elementNumber)
{
    m_elementNumber = elementNumber;
    m_numberOfParticles = m_system->getNumberOfParticles();
    m_numberOfDimensions = m_system->getNumberOfDimensions();
    m_degreesOfFreedom = m_system->getNumberOfFreeDimensions();
    m_numberOfParameters = m_numberOfParticles * m_numberOfParticles;
    m_gradients.zeros(m_system->getMaxParameters());
}

void SimpleJastrow::initializeArrays(const arma::vec positions,
                                     const arma::vec radialVector,
                                     const arma::mat distanceMatrix)
{
    m_positions = positions;
    m_distanceMatrix = distanceMatrix;
    m_probabilityRatio = 1;
    initializePrincipalDistance();
}

void SimpleJastrow::updateArrays(const arma::vec positions,
                                 const arma::vec radialVector,
                                 const arma::mat distanceMatrix,
                                 const arma::uword i)
{
    arma::uword particle = arma::uword(i / m_numberOfDimensions);
    m_positions = positions;
    m_distanceMatrix = distanceMatrix;
    calculateProbabilityRatio(particle);
    updatePrincipalDistance(i, particle);
}

void SimpleJastrow::setArrays()
{
    m_positionsOld = m_positions;
    m_distanceMatrixOld = m_distanceMatrix;
    m_probabilityRatioOld = m_probabilityRatio;
    m_principalDistanceOld = m_principalDistance;
}

void SimpleJastrow::resetArrays()
{
    m_positions = m_positionsOld;
    m_distanceMatrix = m_distanceMatrixOld;
    m_probabilityRatio = m_probabilityRatioOld;
    m_principalDistance = m_principalDistanceOld;
}

void SimpleJastrow::updateParameters(const arma::mat parameters)
{
    arma::vec betaFlatten = parameters.row(m_elementNumber).head(m_numberOfParameters);
    m_beta = arma::reshape(betaFlatten, m_numberOfParticles, m_numberOfParticles);
}

double SimpleJastrow::evaluateRatio()
{
    return m_probabilityRatio;
}

double SimpleJastrow::computeGradient(const arma::uword k)
{
    arma::uword k_p = arma::uword(k / m_numberOfDimensions);
    arma::uword k_d = k % m_numberOfDimensions;

    double derivative = 0;
    for (arma::uword j_p = 0; j_p < k_p; j_p++) {
        arma::uword j = j_p * m_numberOfDimensions + k_d;
        derivative += m_beta(k_p, j_p) * m_principalDistance(k, j);
    }
    for (arma::uword j_p = k_p + 1; j_p < m_numberOfParticles; j_p++) {
        arma::uword j = j_p * m_numberOfDimensions + k_d;
        derivative += m_beta(k_p, j_p) * m_principalDistance(k, j);
    }
    return derivative;
}

double SimpleJastrow::computeLaplacian()
{
    double derivative = 0;
    for (arma::uword k = 0; k < m_degreesOfFreedom; k++) {
        arma::uword k_p = arma::uword(k / m_numberOfDimensions);
        arma::uword k_d = k % m_numberOfDimensions;
        for (arma::uword j_p = 0; j_p < k_p; j_p++) {
            arma::uword j = j_p * m_numberOfDimensions + k_d;
            derivative += m_beta(k_p, j_p)
                          * (1 - m_principalDistance(k, j) * m_principalDistance(k, j))
                          / m_distanceMatrix(k_p, j_p);
        }
        for (arma::uword j_p = k_p + 1; j_p < m_numberOfParticles; j_p++) {
            arma::uword j = j_p * m_numberOfDimensions + k_d;
            derivative += m_beta(k_p, j_p)
                          * (1 - m_principalDistance(k, j) * m_principalDistance(k, j))
                          / m_distanceMatrix(k_p, j_p);
        }
    }
    return derivative;
}

arma::vec SimpleJastrow::computeParameterGradient()
{
    m_gradients.head(m_numberOfParameters) = arma::vectorise(m_distanceMatrix);
    return m_gradients;
}

void SimpleJastrow::initializePrincipalDistance()
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

void SimpleJastrow::updatePrincipalDistance(arma::uword i, arma::uword i_p)
{
    /* Update of the principal distance matrix
     * Arguments:
     * 
     * {arma::uword} i:     The changed coordinate
     * {arma::uword} i_p:   The moved particle
     */
    /*
    arma::uword i_d = i % m_numberOfDimensions;
    for (arma::uword j_p = 0; j_p < i_p; j_p++) {
        arma::uword j = i_d + j_p * m_numberOfDimensions;
        m_principalDistance(i, j) = (m_positions(i) - m_positions(j)) / m_distanceMatrix(i_p, j_p);
        m_principalDistance(j, i) = -m_principalDistance(i, j);
    }
    for (arma::uword j_p = i_p + 1; j_p < m_numberOfParticles; j_p++) {
        arma::uword j = i_d + j_p * m_numberOfDimensions;
        m_principalDistance(i, j) = (m_positions(i) - m_positions(j)) / m_distanceMatrix(i_p, j_p);
        m_principalDistance(j, i) = -m_principalDistance(i, j);
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

void SimpleJastrow::calculateProbabilityRatio(arma::uword i_p)
{
    double ratio = 0;
    for (arma::uword j_p = i_p; j_p < m_numberOfParticles; j_p++) {
        ratio += m_beta(i_p, j_p) * (m_distanceMatrix(i_p, j_p) - m_distanceMatrixOld(i_p, j_p));
    }

    m_probabilityRatio = exp(2 * ratio);
}
