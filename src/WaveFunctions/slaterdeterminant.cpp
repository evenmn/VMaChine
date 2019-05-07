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
    setArrays();

    m_particle  = (unsigned int)(changedCoord/m_numberOfDimensions);
    m_dimension = changedCoord%m_numberOfDimensions;
    m_positions = positions;
    m_positionBlock(m_dimension, m_particle) = m_positions(changedCoord);

    updateSlaterMatrixRow(m_particle);
    for(unsigned int d=(changedCoord-m_dimension); d<(changedCoord+m_numberOfDimensions-m_dimension); d++) {
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

void SlaterDeterminant::resetArrays() {
    m_positions                         = m_positionsOld;
    m_determinantDerivative             = m_determinantDerivativeOld;
    m_determinantSecondDerivative       = m_determinantSecondDerivativeOld;
    m_probabilityRatio                  = m_probabilityRatioOld;
    m_slaterMatrix                      = m_slaterMatrixOld;
    m_slaterMatrixInverse               = m_slaterMatrixInverseOld;
    m_slaterMatrixDer                   = m_slaterMatrixDerOld;
    m_slaterMatrixSecDer                = m_slaterMatrixSecDerOld;
}

void SlaterDeterminant::setArrays() {
    m_positionsOld                      = m_positions;
    m_determinantDerivativeOld          = m_determinantDerivative;
    m_determinantSecondDerivativeOld    = m_determinantSecondDerivative;
    m_probabilityRatioOld               = m_probabilityRatio;
    m_slaterMatrixOld                   = m_slaterMatrix;
    m_slaterMatrixInverseOld            = m_slaterMatrixInverse;
    m_slaterMatrixDerOld                = m_slaterMatrixDer;
    m_slaterMatrixSecDerOld             = m_slaterMatrixSecDer;
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
    setArrays();
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

/*
void SlaterDeterminant::updateSlaterMatrixElement(const int i, const int j) {
    m_slaterMatrix(i,j) = 1;
    for(int k=0; k<m_numberOfDimensions; k++) {
        m_slaterMatrix(i,j) *= m_system->getBasis()->evaluate(m_positions(m_numberOfDimensions*i+k), int(m_listOfStates(j,k)));
    }
}
*/

void SlaterDeterminant::updateSlaterMatrixRow(const unsigned int row) {
    for(unsigned int col=0; col<m_numberOfParticlesHalf; col++) {
        //updateSlaterMatrixElement(row, col);
        m_slaterMatrix(row,col) = m_system->getBasis()->basisElement(col, m_positionBlock.col(row));
    }
}

void SlaterDeterminant::updateSlaterMatrixDerRow(const unsigned int row) {
    unsigned int   particle  = (unsigned int)(row/m_numberOfDimensions);
    unsigned short dimension = row%m_numberOfDimensions;
    for(unsigned int col=0; col<m_numberOfParticlesHalf; col++) {
        m_slaterMatrixDer(row,col) = m_system->getBasis()->basisElementDer(col, dimension, m_positionBlock.col(particle));
    }

    /*
    for(int i=0; i<m_numberOfParticlesHalf; i++) {
        m_slaterMatrixDer(row,i) = m_system->getBasis()->evaluateDerivative(m_positions(row), int(m_listOfStates(i, dimension)));
        for(int j=0; j<m_numberOfDimensions; j++) {
            int m = m_numberOfDimensions * particle + j;
            if(m != row) {
                m_slaterMatrixDer(row,i) *= m_system->getBasis()->evaluate(m_positions(m), int(m_listOfStates(i, j)));
            }
        }
        m_slaterMatrixDer(row,i) = m_system->getBasis()->basisElementDer(i, dimension, m_positionBlock.col(particle));
    }
    */
}

void SlaterDeterminant::updateSlaterMatrixSecDerRow(const unsigned int row) {
    unsigned int   particle  = (unsigned int)(row/m_numberOfDimensions);
    unsigned short dimension = row%m_numberOfDimensions;
    for(unsigned int col=0; col<m_numberOfParticlesHalf; col++) {
        m_slaterMatrixSecDer(row,col) = m_system->getBasis()->basisElementSecDer(col, dimension, m_positionBlock.col(particle));
    }


    /*
    for(int i=0; i<m_numberOfParticlesHalf; i++) {
        m_slaterMatrixSecDer(row,i) = m_system->getBasis()->evaluateSecondDerivative(m_positions(row), int(m_listOfStates(i, dimension)));
        for(int j=0; j<m_numberOfDimensions; j++) {
            int m = m_numberOfDimensions * particle + j;
            if(m != row) {
                m_slaterMatrixSecDer(row,i) *= m_system->getBasis()->evaluate(m_positions(m), int(m_listOfStates(i, j)));
            }
        }
    }
    */
}

double SlaterDeterminant::updateRatio() {
    double R = 0;
    for(unsigned int i=0; i<m_numberOfParticlesHalf; i++) {
        R += m_slaterMatrix(m_particle, i) * m_slaterMatrixInverse(i, m_particle);
    }
    m_probabilityRatio = R*R;
    return R;
}

void SlaterDeterminant::updateSlaterDeterminantDer(unsigned int start, unsigned int end) {
    for(unsigned int i=start*m_numberOfDimensions; i<end*m_numberOfDimensions; i++) {
        m_determinantDerivative(i) = 0;
        m_determinantSecondDerivative(i) = 0;
        unsigned int k = (unsigned int)(i/m_numberOfDimensions);
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

double SlaterDeterminant::evaluateRatio() {
    return m_probabilityRatio;
}

double SlaterDeterminant::computeGradient(const unsigned int k) {
    return m_determinantDerivative(k);
}

double SlaterDeterminant::computeLaplacian() {
    return m_determinantSecondDerivative.sum() - double(m_determinantDerivative.transpose() * m_determinantDerivative);
}

Eigen::VectorXd SlaterDeterminant::computeParameterGradient() {
    return Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);
}
