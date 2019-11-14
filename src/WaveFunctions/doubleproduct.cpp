#include "doubleproduct.h"
#include "../system.h"
#include "wavefunction.h"
#include <cassert>
#include <iostream>

DoubleProduct::DoubleProduct(System *system)
    : WaveFunction(system)
{}

void DoubleProduct::setConstants(const int elementNumber)
{
    m_elementNumber = elementNumber;
    m_numberOfHiddenNodes = m_system->getNumberOfHiddenNodes();
    m_degreesOfFreedom = m_system->getNumberOfFreeDimensions();
    m_numberOfParameters = m_numberOfHiddenNodes * m_degreesOfFreedom + m_numberOfHiddenNodes;
    double sigma = 10 * m_system->getWidth();
    m_sigmaSqrd = sigma * sigma;
    m_sigmaQuad = m_sigmaSqrd * m_sigmaSqrd;
}

void DoubleProduct::updateParameters(Eigen::MatrixXd parameters)
{
    Eigen::VectorXd wFlatten = parameters.row(m_elementNumber)
                                   .segment(m_numberOfHiddenNodes,
                                            m_degreesOfFreedom * m_numberOfHiddenNodes);
    m_W = WaveFunction::reshape(wFlatten, m_degreesOfFreedom, m_numberOfHiddenNodes);

    m_WSqrd = m_W.cwiseAbs2();
    m_b = parameters.row(m_elementNumber).head(m_numberOfHiddenNodes);
}

void DoubleProduct::initializeArrays(const Eigen::VectorXd positions,
                                     const Eigen::VectorXd /*radialVector*/,
                                     const Eigen::MatrixXd /*distanceMatrix*/)
{
    m_positions = positions;
    m_positionBlock = WaveFunction::reshape(positions, m_numberOfParticles, m_numberOfDimensions);
    m_probabilityRatio = 1;

    m_v = Eigen::MatrixXd::Zero(m_numberOfParticles, m_numberOfHiddenNodes);
    m_n = Eigen::MatrixXd::Zero(m_numberOfParticles, m_numberOfHiddenNodes);
    m_p = Eigen::MatrixXd::Zero(m_numberOfParticles, m_numberOfHiddenNodes);
    updateVectors();
}

void DoubleProduct::updateArrays(const Eigen::VectorXd positions,
                                 const Eigen::VectorXd /*radialVector*/,
                                 const Eigen::MatrixXd /*distanceMatrix*/,
                                 const int /*changedCoord*/)
{
    m_positions = positions;
    updateVectors();
    updateRatio();
}

void DoubleProduct::setArrays()
{
    m_positionsOld = m_positions;
    m_vOld = m_v;
    m_nOld = m_n;
    m_pOld = m_p;
    m_probabilityRatioOld = m_probabilityRatio;
}

void DoubleProduct::resetArrays()
{
    m_positions = m_positionsOld;
    m_v = m_vOld;
    m_n = m_nOld;
    m_p = m_pOld;
    m_probabilityRatio = m_probabilityRatioOld;
}

double DoubleProduct::evaluateRatio()
{
    return m_probabilityRatio;
}

double DoubleProduct::computeGradient(const int k)
{
    int particle = int(k % m_numberOfDimensions);
    //double g = double(m_W.row(k) * m_n.row(particle).transpose()) / m_sigmaSqrd;
    return double(m_W.row(k) * m_n.row(particle).transpose()) / m_sigmaSqrd;
}

double DoubleProduct::computeLaplacian()
{
    double result = 0;
    for (int k = 0; k < m_degreesOfFreedom; k++) {
        int particle = int(k % m_numberOfDimensions);
        result += m_WSqrd.row(k) * m_n.row(particle).cwiseProduct(m_p.row(particle)).transpose();
    }
    return result / m_sigmaQuad;
}

Eigen::VectorXd DoubleProduct::computeParameterGradient()
{
    m_gradients = Eigen::VectorXd::Zero(m_system->getMaxParameters());
    Eigen::VectorXd m_nSum = m_n.colwise().sum();
    Eigen::MatrixXd out = m_positions * m_nSum.transpose();
    m_gradients.segment(m_numberOfHiddenNodes, out.size()) = WaveFunction::flatten(out);
    m_gradients.head(m_numberOfHiddenNodes) = m_nSum;
    return m_gradients;
}

void DoubleProduct::updateVectors()
{
    for (int i = 0; i < m_numberOfParticles; i++) {
        m_v.row(i) = m_b
                     + m_positionBlock.row(i)
                           * m_W.block(m_numberOfDimensions * i,
                                       0,
                                       m_numberOfDimensions,
                                       m_numberOfHiddenNodes);
    }
    Eigen::MatrixXd m_e = m_v.array().exp();
    m_p = (m_e + Eigen::MatrixXd::Ones(m_numberOfParticles, m_numberOfHiddenNodes)).cwiseInverse();
    m_n = m_e.cwiseProduct(m_p);
}

void DoubleProduct::updateRatio()
{
    double prod = m_pOld.prod() / m_p.prod();
    m_probabilityRatio = prod * prod;
}
