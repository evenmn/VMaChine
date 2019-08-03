#include "partlyrestricted.h"
#include "../system.h"
#include <cassert>
#include <iostream>

PartlyRestricted::PartlyRestricted(System *system)
    : WaveFunction(system)
{
    m_degreesOfFreedom = m_system->getNumberOfFreeDimensions();
    m_numberOfParameters = m_degreesOfFreedom * m_degreesOfFreedom;
}

void PartlyRestricted::setConstants(const int elementNumber)
{
    m_elementNumber = elementNumber;
    m_gradients = Eigen::VectorXd::Zero(m_system->getMaxParameters());
}

void PartlyRestricted::initializeArrays(const Eigen::VectorXd positions,
                                        const Eigen::VectorXd radialVector,
                                        const Eigen::MatrixXd distanceMatrix)
{
    m_positions = positions;
    m_probabilityRatio = 1;
}

void PartlyRestricted::updateArrays(const Eigen::VectorXd positions,
                                    const Eigen::VectorXd radialVector,
                                    const Eigen::MatrixXd distanceMatrix,
                                    const int changedCoord)
{
    m_positions = positions;
    calculateProbabilityRatio(changedCoord);
}

void PartlyRestricted::setArrays()
{
    m_positionsOld = m_positions;
    m_probabilityRatioOld = m_probabilityRatio;
}

void PartlyRestricted::resetArrays()
{
    m_positions = m_positionsOld;
    m_probabilityRatio = m_probabilityRatioOld;
}

void PartlyRestricted::updateParameters(Eigen::MatrixXd parameters)
{
    m_c = WaveFunction::square(parameters.row(m_elementNumber)
                                    .head(m_numberOfParameters));
}

double PartlyRestricted::evaluateRatio()
{
    return m_probabilityRatio;
}

double PartlyRestricted::computeGradient(const int k)
{
    return 2 * m_c.row(k) * m_positions;
}

double PartlyRestricted::computeLaplacian()
{
    return 2 * m_c.diagonal().sum();
}

Eigen::VectorXd PartlyRestricted::computeParameterGradient()
{
    Eigen::MatrixXd out = m_positions * m_positions.transpose();
    m_gradients.head(out.size()) = WaveFunction::flatten(out);
    return m_gradients;
}

void PartlyRestricted::calculateProbabilityRatio(int i)
{
    double ratio = m_c.row(i) * m_positions;
    m_probabilityRatio = exp(2 * ratio * (m_positions(i) - m_positionsOld(i)));
}
