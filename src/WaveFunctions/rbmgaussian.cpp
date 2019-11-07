#include "rbmgaussian.h"
#include "../system.h"
#include <cassert>
#include <iostream>

RBMGaussian::RBMGaussian(System *system)
    : WaveFunction(system)
{}

void RBMGaussian::setConstants(const int elementNumber)
{
    m_elementNumber = elementNumber;
    m_degreesOfFreedom = m_system->getNumberOfFreeDimensions();
    m_numberOfParameters = m_degreesOfFreedom;
    m_omega = m_system->getFrequency();
    double sigma = m_system->getWidth();
    m_sigmaSqrd = sigma * sigma;
}

void RBMGaussian::updateParameters(const Eigen::MatrixXd parameters)
{
    m_a = (parameters.row(m_elementNumber)).head(m_degreesOfFreedom);
}

void RBMGaussian::initializeArrays(const Eigen::VectorXd positions,
                                   const Eigen::VectorXd /*radialVector*/,
                                   const Eigen::MatrixXd /*distanceMatrix*/)
{
    m_positions = positions;
    m_Xa = positions - m_a;
    m_probabilityRatio = 1;
}

void RBMGaussian::updateArrays(const Eigen::VectorXd positions,
                               const Eigen::VectorXd /*radialVector*/,
                               const Eigen::MatrixXd /*distanceMatrix*/,
                               const int i)
{
    m_positions = positions;
    m_Xa = positions - m_a;
    double expDiff = m_XaOld(i) * m_XaOld(i) - m_Xa(i) * m_Xa(i);
    m_probabilityRatio = exp(expDiff / (2 * m_sigmaSqrd));
}

void RBMGaussian::setArrays()
{
    m_positionsOld = m_positions;
    m_XaOld = m_Xa;
    m_probabilityRatioOld = m_probabilityRatio;
}

void RBMGaussian::resetArrays()
{
    m_positions = m_positionsOld;
    m_Xa = m_XaOld;
    m_probabilityRatio = m_probabilityRatioOld;
}

double RBMGaussian::evaluateRatio()
{
    return m_probabilityRatio;
}

double RBMGaussian::computeGradient(const int k)
{
    return -m_Xa(k) / m_sigmaSqrd;
}

double RBMGaussian::computeLaplacian()
{
    return -m_degreesOfFreedom / m_sigmaSqrd;
}

Eigen::VectorXd RBMGaussian::computeParameterGradient()
{
    m_gradients = Eigen::VectorXd::Zero(m_system->getMaxParameters());
    m_gradients.head(m_degreesOfFreedom) = m_Xa / m_sigmaSqrd;
    return m_gradients;
}
