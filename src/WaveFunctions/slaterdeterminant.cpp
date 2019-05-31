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

void SlaterDeterminant::setConstants(const int elementNumber) {
    m_maxNumberOfParametersPerElement   = m_system->getMaxParameters();
    m_elementNumber                     = elementNumber;
}

void SlaterDeterminant::updateArrays(const Eigen::VectorXd positions, const Eigen::VectorXd radialVector, const Eigen::MatrixXd distanceMatrix, const int changedCoord) {
    m_particle  = int(changedCoord/m_numberOfDimensions);
    m_dimension = changedCoord%m_numberOfDimensions;
    m_positions = positions;
    m_positionBlock(m_dimension, m_particle) = m_positions(changedCoord);
    updateSlaterMatrixRow(m_particle);
    for(int d=(changedCoord-m_dimension); d<(changedCoord+m_numberOfDimensions-m_dimension); d++) {
        updateSlaterMatrixDerRow(d);
        updateSlaterMatrixSecDerRow(d);
    }
    int start = 0, end = m_numberOfParticlesHalf;
    if(m_particle >= m_numberOfParticlesHalf) {
        start   += m_numberOfParticlesHalf;
        end     += m_numberOfParticlesHalf;
    }
    updateSlaterMatrixInverse(start, end);
    updateSlaterDeterminantDerivatives(start, end);
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
    updateSlaterDeterminantDerivatives(0, m_numberOfParticles);
}

void SlaterDeterminant::updateParameters(Eigen::MatrixXd parameters) {}

void SlaterDeterminant::initializeSlaterMatrix() {
    m_slaterMatrix = Eigen::MatrixXd::Ones(m_numberOfParticles, m_numberOfParticlesHalf);
    for(int row=0; row<m_numberOfParticles; row++) {
        updateSlaterMatrixRow(row);
    }
}

void SlaterDeterminant::initializeSlaterMatrixDer() {
    m_slaterMatrixDer = Eigen::MatrixXd::Zero(m_numberOfFreeDimensions, m_numberOfParticlesHalf);
    for(int row=0; row<m_numberOfFreeDimensions; row++) {
        updateSlaterMatrixDerRow(row);
    }
}

void SlaterDeterminant::initializeSlaterMatrixSecDer() {
    m_slaterMatrixSecDer = Eigen::MatrixXd::Zero(m_numberOfFreeDimensions, m_numberOfParticlesHalf);
    for(int row=0; row<m_numberOfFreeDimensions; row++) {
        updateSlaterMatrixSecDerRow(row);
    }
}

void SlaterDeterminant::initializeSlaterMatrixInverse() {
    m_slaterMatrixInverse = Eigen::MatrixXd::Zero(m_numberOfParticlesHalf, m_numberOfParticles);
    m_slaterMatrixInverse.leftCols (m_numberOfParticlesHalf) = m_slaterMatrix.topRows   (m_numberOfParticlesHalf).inverse();
    m_slaterMatrixInverse.rightCols(m_numberOfParticlesHalf) = m_slaterMatrix.bottomRows(m_numberOfParticlesHalf).inverse();
}

void SlaterDeterminant::updateSlaterMatrixRow(const int row) {
    for(int col=0; col<m_numberOfParticlesHalf; col++) {
        m_slaterMatrix(row,col) = m_system->getBasis()->basisElement(col, m_positionBlock.col(row));
    }
}

void SlaterDeterminant::updateSlaterMatrixDerRow(const int row) {
    int particle  = int(row/m_numberOfDimensions);
    int dimension = row%m_numberOfDimensions;
    for(int col=0; col<m_numberOfParticlesHalf; col++) {
        m_slaterMatrixDer(row,col) = m_system->getBasis()->basisElementDer(col, dimension, m_positionBlock.col(particle));
    }
}

void SlaterDeterminant::updateSlaterMatrixSecDerRow(const int row) {
    int particle  = int(row/m_numberOfDimensions);
    int dimension = row%m_numberOfDimensions;
    for(int col=0; col<m_numberOfParticlesHalf; col++) {
        m_slaterMatrixSecDer(row,col) = m_system->getBasis()->basisElementSecDer(col, dimension, m_positionBlock.col(particle));
    }
}

double SlaterDeterminant::updateRatio() {
    double R = 0;
    for(int j=0; j<m_numberOfParticlesHalf; j++) {
        R += m_slaterMatrix(m_particle, j) * m_slaterMatrixInverse(j, m_particle);
    }
    m_probabilityRatio = R*R;
    return R;
}

void SlaterDeterminant::updateSlaterDeterminantDerivatives(int start, int end) {
    for(int i=start*m_numberOfDimensions; i<end*m_numberOfDimensions; i++) {
        m_determinantDerivative(i) = 0;
        m_determinantSecondDerivative(i) = 0;
        int k = int(i/m_numberOfDimensions);
        for(int j=0; j<m_numberOfParticlesHalf; j++) {
            m_determinantDerivative(i)       += m_slaterMatrixDer(i,j)    * m_slaterMatrixInverse(j,k);
            m_determinantSecondDerivative(i) += m_slaterMatrixSecDer(i,j) * m_slaterMatrixInverse(j,k);
        }
    }
}

void SlaterDeterminant::updateSlaterMatrixInverse(int start, int end) {
    double R = updateRatio();
    for(int j=start; j<end; j++) {
        double S = 0;
        for(int l=0; l<m_numberOfParticlesHalf; l++) {
            S += m_slaterMatrix(m_particle, l) * m_slaterMatrixInverse(l,j);
        }
        if(j != m_particle) {
            m_slaterMatrixInverse.col(j) -= S * m_slaterMatrixInverse.col(m_particle) / R;
        }
    }
    m_slaterMatrixInverse.col(m_particle) /= R;
}

double SlaterDeterminant::evaluateRatio() {
    return m_probabilityRatio;
}

double SlaterDeterminant::computeGradient(const int k) {
    return m_determinantDerivative(k);
}

double SlaterDeterminant::computeLaplacian() {
    return m_determinantSecondDerivative.sum() - double(m_determinantDerivative.transpose() * m_determinantDerivative);
}

Eigen::VectorXd SlaterDeterminant::computeParameterGradient() {
    return Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);
}
