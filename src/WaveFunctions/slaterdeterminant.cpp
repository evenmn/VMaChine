#include "slaterdeterminant.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../Basis/basis.h"

SlaterDeterminant::SlaterDeterminant(System* system) :
        WaveFunction(system) {
    m_numberOfFreeDimensions    = m_system->getNumberOfFreeDimensions();
    m_numberOfParticles         = m_system->getNumberOfParticles();
    m_numberOfDimensions        = m_system->getNumberOfDimensions();
    m_totalSpin                 = m_system->getTotalSpin();
    m_numberOfParticlesHalf     = m_numberOfParticles/2.;
    m_numberOfSpinUp            = spinUp();
    m_numberOfSpinDn            = spinDn();
    m_numberOfFreeDimUp         = m_numberOfSpinUp * m_numberOfDimensions;
    m_numberOfFreeDimDn         = m_numberOfSpinDn * m_numberOfDimensions;
}

int SlaterDeterminant::spinUp() {
    return int(m_numberOfParticlesHalf + m_totalSpin);
}

int SlaterDeterminant::spinDn() {
    return int(m_numberOfParticlesHalf - m_totalSpin);
}

void SlaterDeterminant::setConstants(const int elementNumber) {
    m_maxParameters = m_system->getMaxParameters();
    m_elementNumber = elementNumber;
}

void SlaterDeterminant::updateArrays(const Eigen::VectorXd positions, const Eigen::VectorXd radialVector, const Eigen::MatrixXd distanceMatrix, const int changedCoord) {
    m_particle  = int(changedCoord/m_numberOfDimensions);
    m_dimension = changedCoord%m_numberOfDimensions;
    m_positions = positions;
    m_positionBlock(m_dimension, m_particle) = m_positions(changedCoord);
    if(m_particle < m_numberOfSpinUp) {
        updateSlaterMatrixRowUp(m_particle);
        for(int d=(changedCoord-m_dimension); d<(changedCoord+m_numberOfDimensions-m_dimension); d++) {
            updateSlaterMatrixDerRowUp(d);
            updateSlaterMatrixSecDerRowUp(d);
        }
        updateSlaterMatrixInverseUp();
        updateSlaterDeterminantDerivativesUp();
    }
    else {
        m_particle -= m_numberOfSpinUp;
        int chg = changedCoord - m_numberOfFreeDimUp;
        updateSlaterMatrixRowDn(m_particle);
        for(int d=(chg-m_dimension); d<(chg+m_numberOfDimensions-m_dimension); d++) {
            updateSlaterMatrixDerRowDn(d);
            updateSlaterMatrixSecDerRowDn(d);
        }
        updateSlaterMatrixInverseDn();
        updateSlaterDeterminantDerivativesDn();
    }
}

void SlaterDeterminant::setArrays() {
    m_positionsOld                      = m_positions;
    m_positionBlockOld                  = m_positionBlock;
    m_determinantDerivativeOld          = m_determinantDerivative;
    m_determinantSecondDerivativeOld    = m_determinantSecondDerivative;
    m_ratioOld                          = m_ratio;
    m_slaterMatrixUpOld                 = m_slaterMatrixUp;
    m_slaterMatrixInverseUpOld          = m_slaterMatrixInverseUp;
    m_slaterMatrixDerUpOld              = m_slaterMatrixDerUp;
    m_slaterMatrixSecDerUpOld           = m_slaterMatrixSecDerUp;
    m_slaterMatrixDnOld                 = m_slaterMatrixDn;
    m_slaterMatrixInverseDnOld          = m_slaterMatrixInverseDn;
    m_slaterMatrixDerDnOld              = m_slaterMatrixDerDn;
    m_slaterMatrixSecDerDnOld           = m_slaterMatrixSecDerDn;
}

void SlaterDeterminant::resetArrays() {
    m_positions                         = m_positionsOld;
    m_positionBlock                     = m_positionBlockOld;
    m_determinantDerivative             = m_determinantDerivativeOld;
    m_determinantSecondDerivative       = m_determinantSecondDerivativeOld;
    m_ratio                             = m_ratioOld;
    m_slaterMatrixUp                    = m_slaterMatrixUpOld;
    m_slaterMatrixInverseUp             = m_slaterMatrixInverseUpOld;
    m_slaterMatrixDerUp                 = m_slaterMatrixDerUpOld;
    m_slaterMatrixSecDerUp              = m_slaterMatrixSecDerUpOld;
    m_slaterMatrixDn                    = m_slaterMatrixDnOld;
    m_slaterMatrixInverseDn             = m_slaterMatrixInverseDnOld;
    m_slaterMatrixDerDn                 = m_slaterMatrixDerDnOld;
    m_slaterMatrixSecDerDn              = m_slaterMatrixSecDerDnOld;
}

void SlaterDeterminant::initializeArrays(const Eigen::VectorXd positions, const Eigen::VectorXd radialVector, const Eigen::MatrixXd distanceMatrix) {
    m_positions                 = positions;
    m_ratio                     = 1;
    Eigen::Map<Eigen::MatrixXd> positionBlock(m_positions.data(), m_numberOfDimensions, m_numberOfParticles);
    m_positionBlock = positionBlock;
    initializeSlaterMatrix();
    initializeSlaterMatrixDer();
    initializeSlaterMatrixSecDer();
    initializeSlaterMatrixInverse();
    m_determinantDerivative         = Eigen::VectorXd::Zero(m_numberOfFreeDimensions);
    m_determinantSecondDerivative   = Eigen::VectorXd::Zero(m_numberOfFreeDimensions);
    updateSlaterDeterminantDerivativesUp();
    updateSlaterDeterminantDerivativesDn();
}

void SlaterDeterminant::updateParameters(Eigen::MatrixXd parameters) {
    m_system->getBasis()->setParameters(parameters.row(m_elementNumber));
}

void SlaterDeterminant::initializeSlaterMatrix() {
    m_slaterMatrixUp = Eigen::MatrixXd::Zero(m_numberOfSpinUp, m_numberOfSpinUp);
    m_slaterMatrixDn = Eigen::MatrixXd::Zero(m_numberOfSpinDn, m_numberOfSpinDn);
    for(int row=0; row<m_numberOfSpinUp; row++) {
        updateSlaterMatrixRowUp(row);
    }
    for(int row=0; row<m_numberOfSpinDn; row++) {
        updateSlaterMatrixRowDn(row);
    }
}

void SlaterDeterminant::initializeSlaterMatrixDer() {
    m_slaterMatrixDerUp = Eigen::MatrixXd::Zero(m_numberOfFreeDimUp, m_numberOfSpinUp);
    m_slaterMatrixDerDn = Eigen::MatrixXd::Zero(m_numberOfFreeDimDn, m_numberOfSpinDn);
    for(int row=0; row<m_numberOfFreeDimUp; row++) {
        updateSlaterMatrixDerRowUp(row);
    }
    for(int row=0; row<m_numberOfFreeDimDn; row++) {
        updateSlaterMatrixDerRowDn(row);
    }
}

void SlaterDeterminant::initializeSlaterMatrixSecDer() {
    m_slaterMatrixSecDerUp = Eigen::MatrixXd::Zero(m_numberOfFreeDimUp, m_numberOfSpinUp);
    m_slaterMatrixSecDerDn = Eigen::MatrixXd::Zero(m_numberOfFreeDimDn, m_numberOfSpinDn);
    for(int row=0; row<m_numberOfFreeDimUp; row++) {
        updateSlaterMatrixSecDerRowUp(row);
    }
    for(int row=0; row<m_numberOfFreeDimDn; row++) {
        updateSlaterMatrixSecDerRowDn(row);
    }
}

void SlaterDeterminant::initializeSlaterMatrixInverse() {
    m_slaterMatrixInverseUp = m_slaterMatrixUp.inverse();
    m_slaterMatrixInverseDn = m_slaterMatrixDn.inverse();
}

void SlaterDeterminant::updateSlaterMatrixRowUp(const int row) {
    for(int col=0; col<m_numberOfSpinUp; col++) {
        m_slaterMatrixUp(row,col) = m_system->getBasis()->basisElement(col, m_positionBlock.col(row));
    }
}

void SlaterDeterminant::updateSlaterMatrixRowDn(const int row) {
    for(int col=0; col<m_numberOfSpinDn; col++) {
        m_slaterMatrixDn(row,col) = m_system->getBasis()->basisElement(col, m_positionBlock.col(row+m_numberOfSpinUp));
    }
}

void SlaterDeterminant::updateSlaterMatrixDerRowUp(const int row) {
    int particle  = int(row/m_numberOfDimensions);
    int dimension = row%m_numberOfDimensions;
    for(int col=0; col<m_numberOfSpinUp; col++) {
        m_slaterMatrixDerUp(row,col) = m_system->getBasis()->basisElementDer(col, dimension, m_positionBlock.col(particle));
    }
}

void SlaterDeterminant::updateSlaterMatrixDerRowDn(const int row) {
    int particle  = int(row/m_numberOfDimensions) + m_numberOfSpinUp;
    int dimension = row%m_numberOfDimensions;
    for(int col=0; col<m_numberOfSpinDn; col++) {
        m_slaterMatrixDerDn(row,col) = m_system->getBasis()->basisElementDer(col, dimension, m_positionBlock.col(particle));
    }
}

void SlaterDeterminant::updateSlaterMatrixSecDerRowUp(const int row) {
    int particle  = int(row/m_numberOfDimensions);
    int dimension = row%m_numberOfDimensions;
    for(int col=0; col<m_numberOfSpinUp; col++) {
        m_slaterMatrixSecDerUp(row,col) = m_system->getBasis()->basisElementSecDer(col, dimension, m_positionBlock.col(particle));
    }
}

void SlaterDeterminant::updateSlaterMatrixSecDerRowDn(const int row) {
    int particle  = int(row/m_numberOfDimensions) + m_numberOfSpinUp;
    int dimension = row%m_numberOfDimensions;
    for(int col=0; col<m_numberOfSpinDn; col++) {
        m_slaterMatrixSecDerDn(row,col) = m_system->getBasis()->basisElementSecDer(col, dimension, m_positionBlock.col(particle));
    }
}

void SlaterDeterminant::updateRatioUp() {
    m_ratio = m_slaterMatrixUp.row(m_particle) * m_slaterMatrixInverseUp.col(m_particle);
}

void SlaterDeterminant::updateRatioDn() {
    m_ratio = m_slaterMatrixDn.row(m_particle) * m_slaterMatrixInverseDn.col(m_particle);
}

void SlaterDeterminant::updateSlaterDeterminantDerivativesUp() {
    for(int i=0; i<m_numberOfFreeDimUp; i++) {
        int particle = int(i/m_numberOfDimensions);
        m_determinantDerivative(i) = m_slaterMatrixDerUp.row(i) * m_slaterMatrixInverseUp.col(particle);
        m_determinantSecondDerivative(i) = m_slaterMatrixSecDerUp.row(i) * m_slaterMatrixInverseUp.col(particle);
    }
}

void SlaterDeterminant::updateSlaterDeterminantDerivativesDn() {
    for(int i=m_numberOfFreeDimUp; i<m_numberOfFreeDimensions; i++) {
        int particle = int(i/m_numberOfDimensions) - m_numberOfSpinUp;
        m_determinantDerivative(i) = m_slaterMatrixDerDn.row(i-m_numberOfFreeDimUp) * m_slaterMatrixInverseDn.col(particle);
        m_determinantSecondDerivative(i) = m_slaterMatrixSecDerDn.row(i-m_numberOfFreeDimUp) * m_slaterMatrixInverseDn.col(particle);
    }
}

void SlaterDeterminant::updateSlaterMatrixInverseUp() {
    updateRatioUp();
    for(int j=0; j<m_numberOfSpinUp; j++) {
        if(j != m_particle) {
            double S = m_slaterMatrixUp.row(m_particle) * m_slaterMatrixInverseUp.col(j);
            m_slaterMatrixInverseUp.col(j) -= S * m_slaterMatrixInverseUp.col(m_particle) / m_ratio;
        }
    }
    m_slaterMatrixInverseUp.col(m_particle) /= m_ratio;
}

void SlaterDeterminant::updateSlaterMatrixInverseDn() {
    updateRatioDn();
    for(int j=0; j<m_numberOfSpinDn; j++) {
        if(j != m_particle) {
            double S = m_slaterMatrixDn.row(m_particle) * m_slaterMatrixInverseDn.col(j);
            m_slaterMatrixInverseDn.col(j) -= S * m_slaterMatrixInverseDn.col(m_particle) / m_ratio;
        }
    }
    m_slaterMatrixInverseDn.col(m_particle) /= m_ratio;
}

double SlaterDeterminant::evaluateRatio() {
    return m_ratio * m_ratio;
}

double SlaterDeterminant::computeGradient(const int k) {
    return m_determinantDerivative(k);
}

double SlaterDeterminant::computeLaplacian() {
    return m_determinantSecondDerivative.sum() - double(m_determinantDerivative.transpose() * m_determinantDerivative);
}

Eigen::VectorXd SlaterDeterminant::computeParameterGradient() {
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxParameters);
    gradients(0) = m_system->getBasis()->basisElementPar(1,m_positionBlock);
    return gradients;
}
