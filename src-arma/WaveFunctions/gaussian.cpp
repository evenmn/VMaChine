#include "gaussian.h"
#include "../system.h"
#include <cassert>
#include <iostream>

Gaussian::Gaussian(System *system)
    : WaveFunction(system)
{}

void Gaussian::setConstants(const arma::uword elementNumber)
{
    m_elementNumber = elementNumber;
    m_degreesOfFreedom = m_system->getNumberOfFreeDimensions();
    m_omega = m_system->getFrequency();
    m_gradients.zeros(m_system->getMaxParameters());
}

void Gaussian::initializeArrays(const arma::vec positions,
                                const arma::vec radialVector,
                                const arma::mat distanceMatrix)
{
    m_positions = positions;
    m_probabilityRatio = 1;
}

void Gaussian::updateProbabilityRatio(arma::uword changedCoord)
{
    m_probabilityRatio = exp(m_omegalpha
                             * (m_positionsOld(changedCoord) * m_positionsOld(changedCoord)
                                - m_positions(changedCoord) * m_positions(changedCoord)));
}

void Gaussian::updateArrays(const arma::vec positions,
                            const arma::vec radialVector,
                            const arma::mat distanceMatrix,
                            const arma::uword changedCoord)
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

void Gaussian::updateParameters(const arma::mat parameters)
{
    m_alpha = parameters(m_elementNumber, 0);
    m_omegalpha = m_omega * m_alpha;
}

double Gaussian::evaluateRatio()
{
    return m_probabilityRatio;
}

double Gaussian::computeGradient(const arma::uword k)
{
    return -m_omegalpha * m_positions(k);
}

double Gaussian::computeLaplacian()
{
    ;
    return -m_omegalpha * m_degreesOfFreedom;
}

arma::vec Gaussian::computeParameterGradient()
{
    m_gradients(0) = -0.5 * m_omega * arma::sum(arma::square(m_positions));
    return m_gradients;
}
