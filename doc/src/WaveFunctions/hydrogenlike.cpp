#include "hydrogenlike.h"
#include <cassert>
#include "wavefunction.h"
#include "../system.h"

HydrogenLike::HydrogenLike(System* system) :
        WaveFunction(system) {
    m_numberOfParticles = m_system->getNumberOfParticles();
    m_numberOfDimensions = m_system->getNumberOfDimensions();
    m_maxNumberOfParametersPerElement = m_system->getMaxNumberOfParametersPerElement();
}

double HydrogenLike::calculateRadialVectorElement(const Eigen::VectorXd positions, const int par) {

    double sqrtElementWise = 0;
    for(int d=0; d<m_numberOfDimensions; d++) {
        sqrtElementWise += positions(par*m_numberOfDimensions + d) * positions(par*m_numberOfDimensions + d);
    }
    return sqrt(sqrtElementWise);
}

Eigen::VectorXd HydrogenLike::calculateRadialVector(const Eigen::VectorXd positions) {
    Eigen::VectorXd radialVector = Eigen::VectorXd::Zero(m_numberOfParticles);
    for(int i=0; i<m_numberOfParticles; i++) {
        radialVector(i) = calculateRadialVectorElement(positions, i);
    }
    return radialVector;
}

void HydrogenLike::initializeArrays(const Eigen::VectorXd positions) {
    m_positions       = positions;
    m_radialVector    = calculateRadialVector(positions);
    m_ratio = 1;
}

void HydrogenLike::updateArrays(const Eigen::VectorXd positions, const int pRand) {
    m_oldPositions         = m_positions;
    m_positions            = positions;

    m_oldRadialVector      = m_radialVector;
    int particle = int(pRand/m_numberOfDimensions);
    m_radialVector(particle)  = calculateRadialVectorElement(positions, particle);

    m_oldRatio = m_ratio;
    m_ratio = exp( 2 * m_alpha * (m_oldRadialVector(pRand) - m_radialVector(pRand)));
}

void HydrogenLike::resetArrays() {
    m_positions       = m_oldPositions;
    m_radialVector    = m_oldRadialVector;
}

void HydrogenLike::updateParameters(const Eigen::MatrixXd parameters, const int elementNumber) {
    m_elementNumber = elementNumber;
    m_alpha = parameters(m_elementNumber, 0);
}

double HydrogenLike::evaluateRatio() {
    return m_ratio;
}

double HydrogenLike::computeFirstDerivative(const int k) {
    return - m_alpha * m_numberOfParticles;
}

double HydrogenLike::computeSecondDerivative() {;
    return 0;
}

Eigen::VectorXd HydrogenLike::computeFirstEnergyDerivative(const int k) {
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);
    gradients(0) = 0.5 * m_numberOfParticles;
    return gradients;
}

Eigen::VectorXd HydrogenLike::computeSecondEnergyDerivative() {
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);
    return gradients;
}
