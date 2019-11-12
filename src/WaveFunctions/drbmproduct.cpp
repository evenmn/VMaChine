#include "drbmproduct.h"
#include "../system.h"
#include "wavefunction.h"
#include <cassert>
#include <iostream>

DRBMProduct::DRBMProduct(System *system, int numberOfLayers)
    : WaveFunction(system)
{
    m_numberOfLayers = numberOfLayers;
}

template<typename T>
int sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}

void DRBMProduct::setConstants(const int elementNumber)
{
    m_elementNumber = elementNumber;
    m_numberOfHiddenNodes = m_system->getNumberOfHiddenNodes();
    m_degreesOfFreedom = m_system->getNumberOfFreeDimensions();
    m_numberOfParameters = m_numberOfHiddenNodes * (1 + m_numberOfLayers * m_degreesOfFreedom);
    double sigma = 1; //m_system->getWidth();
    m_sigmaSqrd = sigma * sigma;
}

void DRBMProduct::updateGradient()
{
    for (int k = 0; k < m_degreesOfFreedom; k++) {
        for (int j = 0; j < m_numberOfHiddenNodes; j++) {
            m_gradientPart(k, j) = 0;
            for (int n = 1; n < m_numberOfLayers + 1; n++) {
                m_gradientPart(k, j) += n * m_positionsPow(n, k)
                                        * m_W((n - 1) * m_degreesOfFreedom + k, j)
                                        / pow(m_sigmaSqrd, 2 * n);
            }
        }
    }
}

void DRBMProduct::updateLaplacian()
{
    for (int k = 0; k < m_degreesOfFreedom; k++) {
        for (int j = 0; j < m_numberOfHiddenNodes; j++) {
            m_laplacianPart(k, j) = 0;
            for (int n = 0; n < m_numberOfLayers; n++) {
                m_gradientPart(k, j) += n * (n + 1) * m_positionsPow(n, k)
                                        * m_W(n * m_degreesOfFreedom + k, j)
                                        / pow(m_sigmaSqrd, 2 * (n + 1));
            }
        }
    }
}

void DRBMProduct::updateVectors()
{
    m_v = m_b;
    for (int n = 0; n < m_numberOfLayers; n++) {
        m_v += m_positionsPow.row(n + 2)
               * m_W.block(n * m_degreesOfFreedom, 0, m_degreesOfFreedom, m_numberOfHiddenNodes)
               / pow(m_sigmaSqrd, 2 * (n + 1));
    }
    Eigen::VectorXd m_e = m_v.array().exp();
    m_p = (m_e + Eigen::VectorXd::Ones(m_numberOfHiddenNodes)).cwiseInverse();
    m_n = m_e.cwiseProduct(m_p);
}

void DRBMProduct::updateRatio()
{
    double Prod = 1;
    for (int j = 0; j < m_numberOfHiddenNodes; j++) {
        Prod *= m_pOld(j) / m_p(j);
    }
    m_probabilityRatio = Prod * Prod;
}

void DRBMProduct::updateArrays(const Eigen::VectorXd positions,
                               const Eigen::VectorXd radialVector,
                               const Eigen::MatrixXd distanceMatrix,
                               const int changedCoord)
{
    m_positions = positions;
    for (int n = 1; n < m_numberOfLayers + 1; n++) {
        m_positionsPow(n + 1, changedCoord) = pow(m_positions(changedCoord), n);
    }
    updateVectors();
    updateRatio();
    updateGradient();
    updateLaplacian();
}

void DRBMProduct::setArrays()
{
    m_positionsOld = m_positions;
    m_positionsPowOld = m_positionsPow;
    m_gradientPartOld = m_gradientPart;
    m_laplacianPartOld = m_laplacianPart;
    m_vOld = m_v;
    m_nOld = m_n;
    m_pOld = m_p;
    m_probabilityRatioOld = m_probabilityRatio;
}

void DRBMProduct::resetArrays()
{
    m_positions = m_positionsOld;
    m_positionsPow = m_positionsPowOld;
    m_gradientPart = m_gradientPartOld;
    m_laplacianPart = m_laplacianPartOld;
    m_v = m_vOld;
    m_n = m_nOld;
    m_p = m_pOld;
    m_probabilityRatio = m_probabilityRatioOld;
}

void DRBMProduct::initializeArrays(const Eigen::VectorXd positions,
                                   const Eigen::VectorXd radialVector,
                                   const Eigen::MatrixXd distanceMatrix)
{
    m_positions = positions;
    m_positionsPow = Eigen::MatrixXd::Zero(m_numberOfLayers + 2, m_degreesOfFreedom);
    for (int n = 0; n < m_numberOfLayers + 1; n++) {
        m_positionsPow.row(n + 1) = m_positions.array().pow(n);
    }
    m_probabilityRatio = 1;

    m_n = Eigen::VectorXd::Zero(m_numberOfHiddenNodes);
    m_p = Eigen::VectorXd::Zero(m_numberOfHiddenNodes);
    m_gradientPart = Eigen::MatrixXd::Zero(m_degreesOfFreedom, m_numberOfHiddenNodes);
    m_laplacianPart = Eigen::MatrixXd::Zero(m_degreesOfFreedom, m_numberOfHiddenNodes);

    updateVectors();
    updateGradient();
    updateLaplacian();
}

void DRBMProduct::updateParameters(Eigen::MatrixXd parameters)
{
    m_b = parameters.row(m_elementNumber).head(m_numberOfHiddenNodes);
    m_W = Eigen::MatrixXd::Zero(m_numberOfLayers * m_degreesOfFreedom, m_numberOfHiddenNodes);
    for (int n = 0; n < m_numberOfLayers; n++) {
        Eigen::VectorXd wFlatten = parameters.row(m_elementNumber)
                                       .segment(m_numberOfHiddenNodes * (1 + n * m_degreesOfFreedom),
                                                m_degreesOfFreedom * m_numberOfHiddenNodes);
        Eigen::Map<Eigen::MatrixXd> W(wFlatten.data(), m_degreesOfFreedom, m_numberOfHiddenNodes);
        m_W.block(n * m_degreesOfFreedom, 0, m_degreesOfFreedom, m_numberOfHiddenNodes) = W;
        //m_W.block(n*m_degreesOfFreedom, 0, m_degreesOfFreedom, m_numberOfHiddenNodes) = WaveFunction::reshape(wFlatten, m_degreesOfFreedom, m_numberOfHiddenNodes);
    }
    //m_W.block(m_degreesOfFreedom, 0, m_degreesOfFreedom, m_numberOfHiddenNodes) = Eigen::MatrixXd::Zero(m_degreesOfFreedom, m_numberOfHiddenNodes);
    //m_W.block(2*m_degreesOfFreedom, 0, m_degreesOfFreedom, m_numberOfHiddenNodes) = Eigen::MatrixXd::Zero(m_degreesOfFreedom, m_numberOfHiddenNodes);
    //std::cout << m_W << std::endl;
}

double DRBMProduct::evaluateRatio()
{
    return m_probabilityRatio;
}

double DRBMProduct::computeGradient(const int k)
{
    return m_gradientPart.row(k) * m_n;
}

double DRBMProduct::computeLaplacian()
{
    double sum = 0;
    for (int k = 0; k < m_degreesOfFreedom; k++) {
        for (int j = 0; j < m_numberOfHiddenNodes; j++) {
            sum += m_n(j)
                   * (m_laplacianPart(k, j) + m_p(j) * m_gradientPart(k, j) * m_gradientPart(k, j));
        }
    }
    return sum;
}

Eigen::VectorXd DRBMProduct::computeParameterGradient()
{
    m_gradients = Eigen::VectorXd::Zero(m_system->getMaxParameters());
    for (int l = 0; l < m_numberOfHiddenNodes; l++) {
        m_gradients(l) = m_n(l);
        for (int m = 0; m < m_degreesOfFreedom; m++) {
            for (int n = 0; n < m_numberOfLayers; n++) {
                int o = l * m_degreesOfFreedom + m
                        + m_numberOfHiddenNodes * (1 + n * m_degreesOfFreedom);
                m_gradients(o) = m_positionsPow(n + 2, m) * m_n(l) / pow(m_sigmaSqrd, 2 * (n + 1));
            }
        }
    }

    /*
    gradients.head(m_numberOfHiddenNodes) = m_n;
    for(int n=0; n<m_numberOfLayers; n++) {
        Eigen::MatrixXd out = m_positionsPow.row(n+2) * m_n.transpose();
        gradients.segment(m_numberOfHiddenNodes, m_numberOfHiddenNodes*m_degreesOfFreedom) = WaveFunction::flatten(out) / pow(m_sigmaSqrd, 2*(n+1));
    }
    */
    return m_gradients;
}
