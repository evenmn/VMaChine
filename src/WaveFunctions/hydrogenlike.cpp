#include "hydrogenlike.h"
#include <cassert>
#include "wavefunction.h"
#include "../system.h"

HydrogenLike::HydrogenLike(System* system) :
        WaveFunction(system) {
    m_numberOfParticles                 = m_system->getNumberOfParticles();
    m_numberOfDimensions                = m_system->getNumberOfDimensions();
    m_numberOfFreeDimensions            = m_system->getNumberOfFreeDimensions();
    m_Z                                 = m_system->getAtomicNumber();
}

void HydrogenLike::initializeArrays(const Eigen::VectorXd positions, const Eigen::VectorXd radialVector, const Eigen::MatrixXd distanceMatrix) {
    m_positions         = positions;
    m_radialVector      = radialVector;
    m_probabilityRatio  = 1;

    setArrays();
}

void HydrogenLike::updateArrays(const Eigen::VectorXd positions, const Eigen::VectorXd radialVector, const Eigen::MatrixXd distanceMatrix, const int changedCoord) {
    int particle = int(changedCoord/m_numberOfDimensions);

    setArrays();

    m_positions                 = positions;
    m_radialVector              = radialVector;
    m_probabilityRatio          = exp( 2 * m_Z * m_alpha * (m_radialVectorOld(particle) - m_radialVector(particle)));
}

void HydrogenLike::setArrays() {
    m_positionsOld          = m_positions;
    m_radialVectorOld       = m_radialVector;
    m_probabilityRatioOld   = m_probabilityRatio;
}

void HydrogenLike::resetArrays() {
    m_positions             = m_positionsOld;
    m_radialVector          = m_radialVectorOld;
    m_probabilityRatio      = m_probabilityRatioOld;
}

void HydrogenLike::updateParameters(const Eigen::MatrixXd parameters, const int elementNumber) {
    m_elementNumber                     = elementNumber;
    m_maxNumberOfParametersPerElement   = m_system->getMaxNumberOfParametersPerElement();
    m_alpha                             = parameters(m_elementNumber, 0);
}

double HydrogenLike::evaluateRatio() {
    return m_probabilityRatio;
}

double HydrogenLike::computeGradient(const int k) {
    if(k%m_numberOfDimensions == 0) {
        return - m_alpha * m_Z;
    }
    else {
        return 0;
    }
}

double HydrogenLike::computeLaplacian() {
    double derivative = m_radialVector.cwiseInverse().sum();
    return - 2 * m_alpha * m_Z * derivative;
}

Eigen::VectorXd HydrogenLike::computeParameterGradient() {
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);
    gradients(0) = - m_Z * m_radialVector.sum();
    return gradients;
}
