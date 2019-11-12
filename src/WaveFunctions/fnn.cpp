#include "fnn.h"
#include "../system.h"
#include <cassert>
#include <iostream>

FNN::FNN(System *system, int numberOfHiddenUnits, double (*f1)(double), double (*f2)(double))
    : WaveFunction(system)
{
    m_numberOfHiddenUnits = numberOfHiddenUnits;
    m_f1 = f1;
    m_f2 = f2;
}

void FNN::setConstants(const int elementNumber)
{
    m_elementNumber = elementNumber;
    m_degreesOfFreedom = m_system->getNumberOfFreeDimensions();
    m_omega = m_system->getFrequency();
}

void FNN::initializeArrays(const Eigen::VectorXd positions,
                           const Eigen::VectorXd /*radialVector*/,
                           const Eigen::MatrixXd /*distanceMatrix*/)
{
    m_positions = positions;
    m_probabilityRatio = 1;
}

void FNN::updateProbabilityRatio(int changedCoord)
{
    double psiNew = m_probabilityRatio = exp(
        m_omegalpha
        * (m_positionsOld(changedCoord) * m_positionsOld(changedCoord)
           - m_positions(changedCoord) * m_positions(changedCoord)));
}

void FNN::updateArrays(const Eigen::VectorXd positions,
                       const Eigen::VectorXd /*radialVector*/,
                       const Eigen::MatrixXd /*distanceMatrix*/,
                       const int changedCoord)
{
    m_positions = positions;
    updateProbabilityRatio(changedCoord);
}

void FNN::setArrays()
{
    m_positionsOld = m_positions;
    m_probabilityRatioOld = m_probabilityRatio;
}

void FNN::resetArrays()
{
    m_positions = m_positionsOld;
    m_probabilityRatio = m_probabilityRatioOld;
}

void FNN::updateParameters(const Eigen::MatrixXd parameters)
{
    m_alpha = parameters(m_elementNumber, 0);
    m_omegalpha = m_omega * m_alpha;
}

double FNN::evaluateRatio()
{
    return m_probabilityRatio;
}

double FNN::computeGradient(const int k)
{
    return -m_omegalpha * m_positions(k);
}

double FNN::computeLaplacian()
{
    ;
    return -m_omegalpha * m_degreesOfFreedom;
}

Eigen::VectorXd FNN::computeParameterGradient()
{
    m_gradients = Eigen::VectorXd::Zero(m_system->getMaxParameters());
    m_gradients(0) = -0.5 * m_omega * m_positions.cwiseAbs2().sum();
    return m_gradients;
}
