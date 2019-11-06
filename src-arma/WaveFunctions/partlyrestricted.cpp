#include "partlyrestricted.h"
#include "../system.h"
#include <cassert>
#include <iostream>

PartlyRestricted::PartlyRestricted(System *system)
    : WaveFunction(system)
{}

void PartlyRestricted::setConstants(const arma::uword elementNumber)
{
    m_elementNumber = elementNumber;
    m_degreesOfFreedom = m_system->getNumberOfFreeDimensions();
    m_numberOfParameters = m_degreesOfFreedom * m_degreesOfFreedom;
    m_gradients.zeros(m_system->getMaxParameters());
}

void PartlyRestricted::initializeArrays(const arma::vec positions,
                                        const arma::vec radialVector,
                                        const arma::mat distanceMatrix)
{
    m_positions = positions;
    m_probabilityRatio = 1;
}

void PartlyRestricted::updateArrays(const arma::vec positions,
                                    const arma::vec radialVector,
                                    const arma::mat distanceMatrix,
                                    const arma::uword changedCoord)
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

void PartlyRestricted::updateParameters(arma::mat parameters)
{
    m_c = arma::reshape(parameters.row(m_elementNumber)
                                    .head(m_numberOfParameters), m_numberOfParticles, m_numberOfParticles);
}

double PartlyRestricted::evaluateRatio()
{
    return m_probabilityRatio;
}

double PartlyRestricted::computeGradient(const arma::uword k)
{
    return 2 * arma::dot(m_c.row(k), m_positions);
}

double PartlyRestricted::computeLaplacian()
{
    return 2 * arma::trace(m_c);
}

arma::vec PartlyRestricted::computeParameterGradient()
{
    arma::mat out = m_positions * arma::trans(m_positions);
    m_gradients.head(out.size()) = arma::vectorise(out);
    return m_gradients;
}

void PartlyRestricted::calculateProbabilityRatio(arma::uword i)
{
    double ratio = arma::dot(m_c.row(i), m_positions);
    m_probabilityRatio = exp(2 * ratio * (m_positions(i) - m_positionsOld(i)));
}
