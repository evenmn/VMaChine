#include "gaussianspherical.h"
#include "../system.h"
#include <cassert>
#include <iostream>

GaussianSpherical::GaussianSpherical(System *system)
    : WaveFunction(system)
{
    m_degreesOfFreedom = m_system->getNumberOfFreeDimensions();
    m_omega = m_system->getFrequency();
}

void GaussianSpherical::setConstants(const int elementNumber)
{
    m_maxParameters = m_system->getMaxParameters();
    m_elementNumber = elementNumber;
}

void GaussianSpherical::initializeArrays(const Eigen::VectorXd positions,
                                         const Eigen::VectorXd radialVector,
                                         const Eigen::MatrixXd distanceMatrix)
{
    m_positions = positions;
    m_radialVector = radialVector;
    m_probabilityRatio = 1;
}

void GaussianSpherical::updateProbabilityRatio(int changedCoord)
{
    m_probabilityRatio = exp(m_omega * m_alpha
                             * (m_radialVectorOld(changedCoord) * m_radialVectorOld(changedCoord)
                                - m_radialVector(changedCoord) * m_radialVector(changedCoord)));
}

void GaussianSpherical::updateArrays(const Eigen::VectorXd positions,
                                     const Eigen::VectorXd radialVector,
                                     const Eigen::MatrixXd distanceMatrix,
                                     const int changedCoord)
{
    m_positions = positions;
    m_radialVector = radialVector;
    updateProbabilityRatio(changedCoord);
}

void GaussianSpherical::setArrays()
{
    m_positionsOld = m_positions;
    m_probabilityRatioOld = m_probabilityRatio;
    m_radialVectorOld = m_radialVector;
}

void GaussianSpherical::resetArrays()
{
    m_positions = m_positionsOld;
    m_probabilityRatio = m_probabilityRatioOld;
    m_radialVector = m_radialVectorOld;
}

void GaussianSpherical::updateParameters(const Eigen::MatrixXd parameters)
{
    m_alpha = parameters(m_elementNumber, 0);
}

double GaussianSpherical::evaluateRatio()
{
    return m_probabilityRatio;
}

double GaussianSpherical::computeGradient(const int k)
{
    return -m_omega * m_alpha * m_positions(k);
}

double GaussianSpherical::computeLaplacian()
{
    ;
    return -m_omega * m_alpha * m_degreesOfFreedom;
}

Eigen::VectorXd GaussianSpherical::computeParameterGradient()
{
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxParameters);
    gradients(0) = -0.5 * m_omega * m_positions.cwiseAbs2().sum();
    return gradients;
}
