#include "slaterdeterminant.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../Basis/basis.h"

SlaterDeterminant::SlaterDeterminant(System* system) :
        WaveFunction(system) {
    m_numberOfFreeDimensions            = m_system->getNumberOfFreeDimensions();
    m_numberOfParticles                 = m_system->getNumberOfParticles();
    m_numberOfDimensions                = m_system->getNumberOfDimensions();
    m_numberOfParticlesHalf             = m_numberOfParticles/2;
}

void SlaterDeterminant::updateArrays(const Eigen::VectorXd positions, const Eigen::VectorXd radialVector, const Eigen::MatrixXd distanceMatrix, const unsigned int changedCoord) {
    m_particle  = unsigned(changedCoord/m_numberOfDimensions);
    m_dimension = changedCoord%m_numberOfDimensions;
    m_positions = positions;
    m_particleHalf = m_particle/2;
    m_positionBlock(m_dimension, m_particle) = m_positions(changedCoord);
    updateSlaterMatrixRow(m_particle);
    for(unsigned int d=m_particleHalf; d<m_particleHalf+m_numberOfDimensions; d++) {
        updateSlaterMatrixDerRow(d);
        updateSlaterMatrixSecDerRow(d);
    }
    unsigned int start = 0, end = m_numberOfParticlesHalf;
    if(m_particle >= m_numberOfParticlesHalf) {
        start   += m_numberOfParticlesHalf;
        end     += m_numberOfParticlesHalf;
    }
    updateSlaterMatrixInverse (start, end);
    updateSlaterDeterminantDer(start, end);
}

void SlaterDeterminant::setArrays() {
    m_positionsOld                      = m_positions;
    m_positionBlockOld                  = m_positionBlock;
    m_determinantDerivativeOld          = m_determinantDerivative;
    m_determinantSecondDerivativeOld    = m_determinantSecondDerivative;
    m_probabilityRatioOld               = m_probabilityRatio;
    m_slaterMatrixOld                   = m_slaterMatrix;
    m_slaterMatrixInverseOld            = m_slaterMatrixInverse;
    m_slaterMatrixDerOld                = m_slaterMatrixDer;
    m_slaterMatrixSecDerOld             = m_slaterMatrixSecDer;
}

void SlaterDeterminant::resetArrays() {
    m_positions                         = m_positionsOld;
    m_positionBlock                     = m_positionBlockOld;
    m_determinantDerivative             = m_determinantDerivativeOld;
    m_determinantSecondDerivative       = m_determinantSecondDerivativeOld;
    m_probabilityRatio                  = m_probabilityRatioOld;
    m_slaterMatrix                      = m_slaterMatrixOld;
    m_slaterMatrixInverse               = m_slaterMatrixInverseOld;
    m_slaterMatrixDer                   = m_slaterMatrixDerOld;
    m_slaterMatrixSecDer                = m_slaterMatrixSecDerOld;
}

void SlaterDeterminant::initializeArrays(const Eigen::VectorXd positions, const Eigen::VectorXd radialVector, const Eigen::MatrixXd distanceMatrix) {
    // Set matrices
    m_positions                 = positions;
    m_probabilityRatio          = 1;
    Eigen::Map<Eigen::MatrixXd> positionBlock(m_positions.data(), m_numberOfDimensions, m_numberOfParticles);
    m_positionBlock = positionBlock;
    initializeSlaterMatrix();
    initializeSlaterMatrixDer();
    initializeSlaterMatrixSecDer();
    initializeSlaterMatrixInverse();
    m_determinantDerivative         = Eigen::VectorXd::Zero(m_numberOfFreeDimensions);
    m_determinantSecondDerivative   = Eigen::VectorXd::Zero(m_numberOfFreeDimensions);
    updateSlaterDeterminantDer(0, m_numberOfParticles);
}

void SlaterDeterminant::updateParameters(const Eigen::MatrixXd parameters, const unsigned short elementNumber) {
    m_maxNumberOfParametersPerElement   = m_system->getMaxNumberOfParametersPerElement();
    m_elementNumber                     = elementNumber;
}

void SlaterDeterminant::initializeSlaterMatrix() {
    m_slaterMatrix = Eigen::MatrixXd::Ones(m_numberOfParticles, m_numberOfParticlesHalf);
    for(unsigned int row=0; row<m_numberOfParticles; row++) {
        updateSlaterMatrixRow(row);
    }
}

void SlaterDeterminant::initializeSlaterMatrixDer() {
    m_slaterMatrixDer = Eigen::MatrixXd::Zero(m_numberOfFreeDimensions, m_numberOfParticlesHalf);
    for(unsigned int row=0; row<m_numberOfFreeDimensions; row++) {
        updateSlaterMatrixDerRow(row);
    }
}

void SlaterDeterminant::initializeSlaterMatrixSecDer() {
    m_slaterMatrixSecDer = Eigen::MatrixXd::Zero(m_numberOfFreeDimensions, m_numberOfParticlesHalf);
    for(unsigned int row=0; row<m_numberOfFreeDimensions; row++) {
        updateSlaterMatrixSecDerRow(row);
    }
}

void SlaterDeterminant::initializeSlaterMatrixInverse() {
    m_slaterMatrixInverse = Eigen::MatrixXd::Zero(m_numberOfParticlesHalf, m_numberOfParticles);
    m_slaterMatrixInverse.leftCols (m_numberOfParticlesHalf) = m_slaterMatrix.topRows   (m_numberOfParticlesHalf).inverse();
    m_slaterMatrixInverse.rightCols(m_numberOfParticlesHalf) = m_slaterMatrix.bottomRows(m_numberOfParticlesHalf).inverse();
}

void SlaterDeterminant::updateSlaterMatrixRow(const unsigned int row) {
    for(unsigned int col=0; col<m_numberOfParticlesHalf; col++) {
        m_slaterMatrix(row,col) = m_system->getBasis()->basisElement(col, m_positionBlock.col(row));
    }
}

void SlaterDeterminant::updateSlaterMatrixDerRow(const unsigned int row) {
    unsigned int   particle  = unsigned(row/m_numberOfDimensions);
    unsigned short dimension = row%m_numberOfDimensions;
    for(unsigned int col=0; col<m_numberOfParticlesHalf; col++) {
        m_slaterMatrixDer(row,col) = m_system->getBasis()->basisElementDer(col, dimension, m_positionBlock.col(particle));
    }
}

void SlaterDeterminant::updateSlaterMatrixSecDerRow(const unsigned int row) {
    unsigned int   particle  = unsigned(row/m_numberOfDimensions);
    unsigned short dimension = row%m_numberOfDimensions;
    for(unsigned int col=0; col<m_numberOfParticlesHalf; col++) {
        m_slaterMatrixSecDer(row,col) = m_system->getBasis()->basisElementSecDer(col, dimension, m_positionBlock.col(particle));
    }
}

void SlaterDeterminant::updateSlaterDeterminantDer(unsigned int start, unsigned int end) {
    for(unsigned int i=start*m_numberOfDimensions; i<end*m_numberOfDimensions; i++) {
        m_determinantDerivative(i) = 0;
        m_determinantSecondDerivative(i) = 0;
        unsigned int k = unsigned(i/m_numberOfDimensions);
        for(unsigned int j=0; j<m_numberOfParticlesHalf; j++) {
            m_determinantDerivative(i)       += m_slaterMatrixDer(i,j)    * m_slaterMatrixInverse(j,k);
            m_determinantSecondDerivative(i) += m_slaterMatrixSecDer(i,j) * m_slaterMatrixInverse(j,k);
        }
    }
}

void SlaterDeterminant::updateSlaterMatrixInverse(const unsigned int start, const unsigned int end) {
    double R = updateRatio();
    for(unsigned int i=start; i<end; i++) {
        double S = 0;
        for(unsigned int j=0; j<m_numberOfParticlesHalf; j++) {
            S += m_slaterMatrix(m_particle, j) * m_slaterMatrixInverse(j,i);
        }
        if(i != m_particle) {
            m_slaterMatrixInverse.col(i) -= S * m_slaterMatrixInverse.col(m_particle) / R;
        }
    }
    m_slaterMatrixInverse.col(m_particle) /= R;
}

double SlaterDeterminant::updateRatio() {
    double R = 0;
    for(unsigned int i=0; i<m_numberOfParticlesHalf; i++) {
        R += m_slaterMatrix(m_particle, i) * m_slaterMatrixInverse(i, m_particle);
    }
    m_probabilityRatio = R*R;
    return R;
}

double SlaterDeterminant::evaluateRatio() {
    return m_probabilityRatio;
}

double SlaterDeterminant::computeGradient(const unsigned int k) {
    return m_determinantDerivative(k);
}

double SlaterDeterminant::computeLaplacian() {
    //std::cout << m_determinantSecondDerivative.sum() << std::endl;
    return m_determinantSecondDerivative.sum() - m_determinantDerivative.cwiseAbs2().sum();
    //return - m_determinantDerivative.cwiseAbs2().sum();
}

Eigen::VectorXd SlaterDeterminant::computeParameterGradient() {
    return Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);
}
