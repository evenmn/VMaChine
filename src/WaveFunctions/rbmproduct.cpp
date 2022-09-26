#include "rbmproduct.h"
#include "../system.h"
#include "wavefunction.h"

RBMProduct::RBMProduct(System *system)
    : WaveFunction(system)
{}

void RBMProduct::setConstants(const int elementNumber)
{
    m_elementNumber = elementNumber;
    m_numberOfHiddenNodes = m_system->getNumberOfHiddenUnits();
    m_degreesOfFreedom = m_system->getNumberOfFreeDimensions();
    m_numberOfParameters = m_numberOfHiddenNodes * m_degreesOfFreedom + m_numberOfHiddenNodes;
    double sigma = 10 * m_system->getWidth();
    m_sigmaSqrd = sigma * sigma;
    m_sigmaQuad = m_sigmaSqrd * m_sigmaSqrd;
    m_sorting = m_system->getSorting();
}

void RBMProduct::updateParameters(Eigen::MatrixXd parameters)
{
    Eigen::VectorXd wFlatten = parameters.row(m_elementNumber)
                                   .segment(m_numberOfHiddenNodes,
                                            m_degreesOfFreedom * m_numberOfHiddenNodes);
    m_W = WaveFunction::reshape(wFlatten, m_degreesOfFreedom, m_numberOfHiddenNodes);

    m_WSqrd = m_W.cwiseAbs2();
    m_b = parameters.row(m_elementNumber).head(m_numberOfHiddenNodes);
}

void RBMProduct::initializeArrays(const Eigen::VectorXd positions,
                                  const Eigen::VectorXd /*radialVector*/,
                                  const Eigen::MatrixXd /*distanceMatrix*/)
{
    m_positions = positions;
    if (m_sorting) argsort(m_positions, m_sortidx);

    m_probabilityRatio = 1;

    m_n = Eigen::VectorXd::Zero(m_numberOfHiddenNodes);
    m_p = Eigen::VectorXd::Zero(m_numberOfHiddenNodes);
    updateVectors();
}

void RBMProduct::updateArrays(const Eigen::VectorXd positions,
                              const Eigen::VectorXd /*radialVector*/,
                              const Eigen::MatrixXd /*distanceMatrix*/,
                              const int /*changedCoord*/)
{
    m_positions = positions;
    if (m_sorting) argsort(m_positions, m_sortidx);
    updateVectors();
    updateRatio();
}

void RBMProduct::setArrays()
{
    m_positionsOld = m_positions;
    m_vOld = m_v;
    m_nOld = m_n;
    m_pOld = m_p;
    m_probabilityRatioOld = m_probabilityRatio;
}

void RBMProduct::resetArrays()
{
    m_positions = m_positionsOld;
    m_v = m_vOld;
    m_n = m_nOld;
    m_p = m_pOld;
    m_probabilityRatio = m_probabilityRatioOld;
}

double RBMProduct::evaluateRatio()
{
    return m_probabilityRatio;
}

double RBMProduct::computeGradient(const int k)
{
    if (m_sorting) {
        return double(m_W.row(m_sortidx(k)) * m_n) / m_sigmaSqrd;
    }
    else {
        return double(m_W.row(k) * m_n) / m_sigmaSqrd;
    }
}

double RBMProduct::computeLaplacian()
{
    return (m_WSqrd * m_p.cwiseProduct(m_n)).sum() / m_sigmaQuad;
}

Eigen::VectorXd RBMProduct::computeParameterGradient()
{
    m_gradients = Eigen::VectorXd::Zero(m_system->getMaxParameters());
    Eigen::MatrixXd out = m_positions * m_n.transpose();
    m_gradients.segment(m_numberOfHiddenNodes, out.size()) = WaveFunction::flatten(out);
    m_gradients.head(m_numberOfHiddenNodes) = m_n;
    return m_gradients;
}

void RBMProduct::updateVectors()
{
    m_v = m_b + m_W.transpose() * m_positions / m_sigmaSqrd;
    Eigen::VectorXd m_e = m_v.array().exp();
    m_p = (m_e + Eigen::VectorXd::Ones(m_numberOfHiddenNodes)).cwiseInverse();
    m_n = m_e.cwiseProduct(m_p);
}

void RBMProduct::updateRatio()
{
    double prod = m_pOld.prod() / m_p.prod();
    m_probabilityRatio = prod * prod;
}
