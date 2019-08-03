#include "hydrogenlike.h"
#include "../system.h"
#include "wavefunction.h"
#include <cassert>

HydrogenLike::HydrogenLike(System *system)
    : WaveFunction(system)
{
    m_numberOfParticles = m_system->getNumberOfParticles();
    m_numberOfDimensions = m_system->getNumberOfDimensions();
    m_degreesOfFreedom = m_system->getNumberOfFreeDimensions();
    m_Z = m_system->getAtomicNumber();
}

void HydrogenLike::setConstants(const int elementNumber)
{
    m_elementNumber = elementNumber;
    m_gradients = Eigen::VectorXd::Zero(m_system->getMaxParameters());
}

void HydrogenLike::initializeArrays(const Eigen::VectorXd positions,
                                    const Eigen::VectorXd radialVector,
                                    const Eigen::MatrixXd distanceMatrix)
{
    m_positions = positions;
    m_radialVector = radialVector;
    m_probabilityRatio = 1;
}

void HydrogenLike::updateArrays(const Eigen::VectorXd positions,
                                const Eigen::VectorXd radialVector,
                                const Eigen::MatrixXd distanceMatrix,
                                const int changedCoord)
{
    int particle = int(changedCoord / m_numberOfDimensions);

    m_positions = positions;
    m_radialVector = radialVector;
    m_probabilityRatio = exp(2 * m_Z * m_alpha
                             * (m_radialVectorOld(particle) - m_radialVector(particle)));
}

void HydrogenLike::setArrays()
{
    m_positionsOld = m_positions;
    m_radialVectorOld = m_radialVector;
    m_probabilityRatioOld = m_probabilityRatio;
}

void HydrogenLike::resetArrays()
{
    m_positions = m_positionsOld;
    m_radialVector = m_radialVectorOld;
    m_probabilityRatio = m_probabilityRatioOld;
}

void HydrogenLike::updateParameters(const Eigen::MatrixXd parameters)
{
    m_alpha = parameters(m_elementNumber, 0);
}

double HydrogenLike::evaluateRatio()
{
    return m_probabilityRatio;
}

double HydrogenLike::computeGradient(const int k)
{
    return -m_alpha * m_Z * m_positions(k) / m_radialVector(int(k / m_numberOfDimensions));
}

double HydrogenLike::computeLaplacian()
{
    double sum = 0;
    for (int i = 0; i < m_numberOfParticles; i++) {
        sum -= 2 / m_radialVector(i);
    }
    return m_alpha * m_Z * sum;
}

Eigen::VectorXd HydrogenLike::computeParameterGradient()
{
    m_gradients(0) = -m_Z * m_radialVector.sum();
    return m_gradients;
}
