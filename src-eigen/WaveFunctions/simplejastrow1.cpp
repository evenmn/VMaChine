#include "simplejastrow.h"
#include "../system.h"
#include "wavefunction.h"
#include <cassert>
#include <iostream>

SimpleJastrow::SimpleJastrow(System *system)
    : WaveFunction(system)
{
    m_numberOfParticles = m_system->getNumberOfParticles();
    m_numberOfDimensions = m_system->getNumberOfDimensions();
    m_degreesOfFreedom = m_system->getNumberOfFreeDimensions();
    m_numberOfParameters = m_numberOfParticles * m_numberOfParticles;
}

void SimpleJastrow::setConstants(const int elementNumber)
{
    m_elementNumber = elementNumber;
    m_gradients = Eigen::VectorXd::Zero(m_system->getMaxParameters());
}

void SimpleJastrow::initializeArrays(const Eigen::VectorXd positions,
                                     const Eigen::VectorXd radialVector,
                                     const Eigen::MatrixXd distanceMatrix)
{
    m_positions = positions;
    m_distanceMatrix = distanceMatrix;
    m_probabilityRatio = 1;
    initializePrincipalDistance();
}

void SimpleJastrow::updateArrays(const Eigen::VectorXd positions,
                                 const Eigen::VectorXd radialVector,
                                 const Eigen::MatrixXd distanceMatrix,
                                 const int i)
{
    m_particle = int(i / m_numberOfDimensions);
    m_positions = positions;
    m_distanceMatrix = distanceMatrix;
    calculateProbabilityRatio(m_particle);
    updatePrincipalDistance(i);
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

void SimpleJastrow::updateParameters(const Eigen::MatrixXd parameters)
{
    Eigen::VectorXd betaFlatten = parameters.row(m_elementNumber).head(m_numberOfParameters);
    m_beta = WaveFunction::square(betaFlatten);
}

double SimpleJastrow::evaluateRatio()
{
    return m_probabilityRatio;
}

double SimpleJastrow::computeGradient(const int k)
{
    int k_p = int(k / m_numberOfDimensions);
    int k_d = k % m_numberOfDimensions;

    double derivative = 0;
    for (int j_p = 0; j_p < k_p; j_p++) {
        int j = j_p * m_numberOfDimensions + k_d;
        derivative += m_beta(k_p, j_p) * m_principalDistance(k, j);
    }
    for (int j_p = k_p + 1; j_p < m_numberOfParticles; j_p++) {
        int j = j_p * m_numberOfDimensions + k_d;
        derivative += m_beta(k_p, j_p) * m_principalDistance(k, j);
    }
    return derivative;
}

double SimpleJastrow::computeLaplacian()
{
    double derivative = 0;
    for (int k = 0; k < m_degreesOfFreedom; k++) {
        int k_p = int(k / m_numberOfDimensions);
        int k_d = k % m_numberOfDimensions;
        for (int j_p = 0; j_p < k_p; j_p++) {
            int j = j_p * m_numberOfDimensions + k_d;
            derivative += m_beta(k_p, j_p)
                          * (1 - m_principalDistance(k, j) * m_principalDistance(k, j))
                          / m_distanceMatrix(k_p, j_p);
        }
        for (int j_p = k_p + 1; j_p < m_numberOfParticles; j_p++) {
            int j = j_p * m_numberOfDimensions + k_d;
            derivative += m_beta(k_p, j_p)
                          * (1 - m_principalDistance(k, j) * m_principalDistance(k, j))
                          / m_distanceMatrix(k_p, j_p);
        }
    }
    return derivative;
}

Eigen::VectorXd SimpleJastrow::computeParameterGradient()
{
    m_gradients.head(m_numberOfParameters) = WaveFunction::flatten(m_distanceMatrix);
    return m_gradients;
}

void SimpleJastrow::initializePrincipalDistance()
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

void SimpleJastrow::updatePrincipalDistance(int i)
{
    int i_d = i % m_numberOfDimensions;
    for (int j_p = 0; j_p < m_particle; j_p++) {
        int j = i_d + j_p * m_numberOfDimensions;
        m_principalDistance(i, j) = (m_positions(i) - m_positions(j))
                                    / m_distanceMatrix(m_particle, j_p);
        m_principalDistance(j, i) = -m_principalDistance(i, j_p);
    }
    for (int j_p = m_particle + 1; j_p < m_numberOfParticles; j_p++) {
        int j = i_d + j_p * m_numberOfDimensions;
        m_principalDistance(i, j) = (m_positions(i) - m_positions(j))
                                    / m_distanceMatrix(m_particle, j_p);
        m_principalDistance(j, i) = -m_principalDistance(i, j_p);
    }
}

void SimpleJastrow::calculateProbabilityRatio(int i_p)
{
    double ratio = 0;
    for (int j_p = i_p; j_p < m_numberOfParticles; j_p++) {
        ratio += m_beta(i_p, j_p) * (m_distanceMatrix(i_p, j_p) - m_distanceMatrixOld(i_p, j_p));
    }
    m_probabilityRatio = exp(2 * ratio);
}