#include "hardcorejastrow.h"
#include "../system.h"
#include "wavefunction.h"

HardCoreJastrow::HardCoreJastrow(System *system, const double radius)
    : WaveFunction(system)
{
    m_sphereRadius = radius;
}

void HardCoreJastrow::setConstants(const int elementNumber)
{
    m_elementNumber = elementNumber;
    m_numberOfParticles = m_system->getNumberOfParticles();
    m_numberOfDimensions = m_system->getNumberOfDimensions();
    m_degreesOfFreedom = m_system->getNumberOfFreeDimensions();
    m_numberOfParameters = 0;
}

void HardCoreJastrow::initializeArrays(const Eigen::VectorXd positions,
                                     const Eigen::VectorXd /*radialVector*/,
                                     const Eigen::MatrixXd distanceMatrix)
{
    m_positions = positions;
    m_distanceMatrix = distanceMatrix;
    m_distanceMatrixSqrd = distanceMatrix.cwiseAbs2();
    m_probabilityRatio = 1;
    initializePrincipalDistance();
}

void HardCoreJastrow::updateArrays(const Eigen::VectorXd positions,
                                 const Eigen::VectorXd /*radialVector*/,
                                 const Eigen::MatrixXd distanceMatrix,
                                 const int i)
{
    int particle = int(i / m_numberOfDimensions);
    m_positions = positions;
    m_distanceMatrix = distanceMatrix;
    m_distanceMatrixSqrd = distanceMatrix.cwiseAbs2();
    calculateProbabilityRatio(particle);
    updatePrincipalDistance(i, particle);
}

void HardCoreJastrow::setArrays()
{
    m_positionsOld = m_positions;
    m_distanceMatrixOld = m_distanceMatrix;
    m_distanceMatrixSqrdOld = m_distanceMatrixSqrd;
    m_probabilityRatioOld = m_probabilityRatio;
    m_principalDistanceOld = m_principalDistance;
}

void HardCoreJastrow::resetArrays()
{
    m_positions = m_positionsOld;
    m_distanceMatrix = m_distanceMatrixOld;
    m_distanceMatrixSqrd = m_distanceMatrixSqrdOld;
    m_probabilityRatio = m_probabilityRatioOld;
    m_principalDistance = m_principalDistanceOld;
}

void HardCoreJastrow::updateParameters(const Eigen::MatrixXd /*parameters*/) {}

double HardCoreJastrow::evaluateRatio()
{
    return m_probabilityRatio;
}

double HardCoreJastrow::computeGradient(const int k)
{
    int k_p = int(k / m_numberOfDimensions);
    int k_d = k % m_numberOfDimensions;

    double derivative = 0;
    for (int j_p = 0; j_p < k_p; j_p++) {
        int j = j_p * m_numberOfDimensions + k_d;
        derivative += m_sphereRadius * m_principalDistance(j, k) / m_distanceMatrix(j_p, k_p);
    }
    for (int j_p = k_p + 1; j_p < m_numberOfParticles; j_p++) {
        int j = j_p * m_numberOfDimensions + k_d;
        derivative += m_sphereRadius * m_principalDistance(k, j) / m_distanceMatrix(k_p, j_p);
    }
    if (std::isnan(derivative)) {
        std::cout << "computeGradient" << std::endl;
        std::cout << m_distanceMatrix << std::endl;
        exit(0);
    }
    return derivative;
}

double HardCoreJastrow::computeLaplacian()
{
    double derivative = 0;
    for (int k = 0; k < m_degreesOfFreedom; k++) {
        int k_p = int(k / m_numberOfDimensions);
        int k_d = k % m_numberOfDimensions;
        for (int j_p = 0; j_p < k_p; j_p++) {
            int j = j_p * m_numberOfDimensions + k_d;
            derivative -= m_sphereRadius / (m_distanceMatrix(j_p, k_p) * m_distanceMatrixSqrd(j_p, k_p))
                          * (1 + 3 * m_principalDistance(j, k) * m_principalDistance(j, k));
        }
        for (int j_p = k_p + 1; j_p < m_numberOfParticles; j_p++) {
            int j = j_p * m_numberOfDimensions + k_d;
            derivative -= m_sphereRadius / (m_distanceMatrix(k_p, j_p) * m_distanceMatrixSqrd(k_p, j_p))
                          * (1 + 3 * m_principalDistance(k, j) * m_principalDistance(k, j));
        }
    }
    if (std::isnan(derivative)) {
        std::cout << "computeLaplacian" << std::endl;
        exit(0);
    }
    return derivative;
}

Eigen::VectorXd HardCoreJastrow::computeParameterGradient()
{
    return Eigen::VectorXd::Zero(m_system->getMaxParameters());
}

void HardCoreJastrow::initializePrincipalDistance()
{
    m_principalDistance = Eigen::MatrixXd::Zero(m_degreesOfFreedom, m_degreesOfFreedom);
    for (int n = 1; n < m_numberOfParticles; n++) {
        int step = n * m_numberOfDimensions;
        for (int j = 0; j < m_degreesOfFreedom - step; j++) {
            int p = int(j / m_numberOfDimensions);
            m_principalDistance(j, j + step) = (m_positions(j) - m_positions(j + step))
                                               / m_distanceMatrix(p, p + n);
        }
    }
}

void HardCoreJastrow::updatePrincipalDistance(int /*i*/, int /*i_p*/)
{
    /* Update of the principal distance matrix
     * Arguments:
     *
     * {int} i:     The changed coordinate
     * {int} i_p:   The moved particle
     */

    for (int n = 1; n < m_numberOfParticles; n++) {
        int step = n * m_numberOfDimensions;
        for (int j = 0; j < m_degreesOfFreedom - step; j++) {
            int p = int(j / m_numberOfDimensions);
            m_principalDistance(j, j + step) = (m_positions(j) - m_positions(j + step))
                                               / m_distanceMatrix(p, p + n);
        }
    }
}

void HardCoreJastrow::calculateProbabilityRatio(int i_p)
{
    double prod = 1;
    for (int j_p = i_p+1; j_p < m_numberOfParticles; j_p++) {
        if (m_distanceMatrix(i_p, j_p) > m_sphereRadius) {
            double new_ = (1 - m_sphereRadius/m_distanceMatrix(i_p, j_p));
            double old_ = (1 - m_sphereRadius/m_distanceMatrixOld(i_p, j_p));
            prod *= new_/old_;
        }
        else {
            prod = 0;
            break;
        }
    }
    if (std::isnan(prod)){
        std::cout << "calculateProbabilityRatio" << std::endl;
        exit(0);
    }
    m_probabilityRatio = prod * prod;
}
