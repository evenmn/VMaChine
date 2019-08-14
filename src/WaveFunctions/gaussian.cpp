#include "gaussian.h"
#include "../system.h"
#include <cassert>
#include <iostream>

Gaussian::Gaussian(System *system)
    : WaveFunction(system)
{
    m_degreesOfFreedom = m_system->getNumberOfFreeDimensions();
    m_omega = m_system->getFrequency();
}

void Gaussian::setConstants(const int elementNumber)
{
    m_elementNumber = elementNumber;
    m_gradients = Eigen::VectorXd::Zero(m_system->getMaxParameters());
}

void Gaussian::initializeArrays(const Eigen::VectorXd positions,
                                const Eigen::VectorXd radialVector,
                                const Eigen::MatrixXd distanceMatrix)
{
    m_positions = positions;
    m_probabilityRatio = 1;
}

void Gaussian::updateProbabilityRatio(int changedCoord)
{
    m_probabilityRatio = exp(m_omegalpha
                             * (m_positionsOld(changedCoord) * m_positionsOld(changedCoord)
                                - m_positions(changedCoord) * m_positions(changedCoord)));
}

void Gaussian::updateArrays(const Eigen::VectorXd positions,
                            const Eigen::VectorXd radialVector,
                            const Eigen::MatrixXd distanceMatrix,
                            const int changedCoord)
{
    m_positions = positions;
    updateProbabilityRatio(changedCoord);
}

void Gaussian::setArrays()
{
    m_positionsOld = m_positions;
    m_probabilityRatioOld = m_probabilityRatio;
}

void Gaussian::resetArrays()
{
    m_positions = m_positionsOld;
    m_probabilityRatio = m_probabilityRatioOld;
}

void Gaussian::updateParameters(const Eigen::MatrixXd parameters)
{
    m_alpha = parameters(m_elementNumber, 0);
    m_omegalpha = m_omega * m_alpha;
}

double Gaussian::evaluateRatio()
{
    return m_probabilityRatio;
}

double Gaussian::computeGradient(const int k)
{
    return -m_omegalpha * m_positions(k);
}

double Gaussian::computeLaplacian()
{
    ;
    return -m_omegalpha * m_degreesOfFreedom;
}

Eigen::VectorXd Gaussian::computeParameterGradient()
{
    m_gradients(0) = -0.5 * m_omega * m_positions.cwiseAbs2().sum();
    return m_gradients;
}
