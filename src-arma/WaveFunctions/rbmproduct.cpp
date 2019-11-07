#include "rbmproduct.h"
#include "../system.h"
#include "wavefunction.h"
#include <cassert>
#include <iostream>

RBMProduct::RBMProduct(System *system)
    : WaveFunction(system)
{}

void RBMProduct::setConstants(const arma::uword elementNumber)
{
    m_elementNumber = elementNumber;
    m_numberOfHiddenNodes = m_system->getNumberOfHiddenNodes();
    m_degreesOfFreedom = m_system->getNumberOfFreeDimensions();
    m_numberOfParameters = m_numberOfHiddenNodes * m_degreesOfFreedom + m_numberOfHiddenNodes;
    double sigma = 10 * m_system->getWidth();
    m_sigmaSqrd = sigma * sigma;
    m_sigmaQuad = m_sigmaSqrd * m_sigmaSqrd;
    m_gradients.zeros(m_system->getMaxParameters());
}

void RBMProduct::updateParameters(arma::mat parameters)
{
    arma::vec wFlatten = parameters.row(m_elementNumber).subvec(m_numberOfHiddenNodes, (1+m_degreesOfFreedom) * m_numberOfHiddenNodes);
    m_W = arma::reshape(wFlatten, m_degreesOfFreedom, m_numberOfHiddenNodes);

    m_WSqrd = arma::square(m_W);
    m_b = parameters.row(m_elementNumber).head(m_numberOfHiddenNodes);
}

void RBMProduct::initializeArrays(const arma::vec positions,
                                  const arma::vec radialVector,
                                  const arma::mat distanceMatrix)
{
    m_positions = positions;
    m_probabilityRatio = 1;

    m_n.zeros(m_numberOfHiddenNodes);
    m_p.zeros(m_numberOfHiddenNodes);
    updateVectors();
}

void RBMProduct::updateArrays(const arma::vec positions,
                              const arma::vec radialVector,
                              const arma::mat distanceMatrix,
                              const arma::uword changedCoord)
{
    m_positions = positions;
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

double RBMProduct::computeGradient(const arma::uword k)
{
    return arma::dot(m_W.row(k), m_n) / m_sigmaSqrd;
}

double RBMProduct::computeLaplacian()
{
    return arma::sum(arma::dot(m_WSqrd, m_p % m_n)) / m_sigmaQuad;
}

arma::vec RBMProduct::computeParameterGradient()
{
    arma::mat out = m_positions * arma::trans(m_n);
    m_gradients.subvec(m_numberOfHiddenNodes, m_numberOfHiddenNodes + out.size()) = arma::vectorise(out);
    m_gradients.head(m_numberOfHiddenNodes) = m_n;
    return m_gradients;
}

void RBMProduct::updateVectors()
{
    m_v = m_b + arma::dot(m_W, m_positions) / m_sigmaSqrd;
    arma::vec m_e = arma::exp(m_v);
    m_p = arma::pow(m_e + 1, -1);
    //arma::ones<arma::vec>(m_numberOfHiddenNodes)
    m_n = m_e % m_p;
}

void RBMProduct::updateRatio()
{
    double prod = arma::prod(m_pOld) / arma::prod(m_p);
    m_probabilityRatio = prod * prod;
}
