#include "slaterdeterminant.h"
#include "../Basis/basis.h"
#include "../system.h"
#include <cassert>
#include <iostream>

SlaterDeterminant::SlaterDeterminant(System *system)
    : WaveFunction(system)
{}

void SlaterDeterminant::setConstants(const arma::uword elementNumber)
{
    m_elementNumber = elementNumber;
    m_degreesOfFreedom = m_system->getNumberOfFreeDimensions();
    m_numberOfParticles = m_system->getNumberOfParticles();
    m_numberOfDimensions = m_system->getNumberOfDimensions();
    m_numberOfParticlesHalf = m_numberOfParticles / 2;
    m_basis = m_system->getBasis();
    m_gradients.zeros(m_system->getMaxParameters());
}

void SlaterDeterminant::updateArrays(const arma::vec positions,
                                     const arma::vec radialVector,
                                     const arma::mat distanceMatrix,
                                     const arma::uword changedCoord)
{
    m_particle = arma::uword(changedCoord / m_numberOfDimensions);
    m_dimension = changedCoord % m_numberOfDimensions;
    m_positions(m_dimension, m_particle) = positions(changedCoord);
    updateSlaterMatrixRow(m_particle);
    for (arma::uword d = (changedCoord - m_dimension);
         d < (changedCoord + m_numberOfDimensions - m_dimension);
         d++) {
        updateSlaterMatrixDerRow(d);
        updateSlaterMatrixSecDerRow(d);
    }
    arma::uword start = 0, end = m_numberOfParticlesHalf;
    if (m_particle >= m_numberOfParticlesHalf) {
        start += m_numberOfParticlesHalf;
        end += m_numberOfParticlesHalf;
    }
    updateSlaterMatrixInverse(start, end);
    updateSlaterDeterminantDerivatives(start, end);
}

void SlaterDeterminant::setArrays()
{
    m_positionsOld = m_positions;
    m_determinantDerivativeOld = m_determinantDerivative;
    m_determinantSecondDerivativeOld = m_determinantSecondDerivative;
    m_ratioOld = m_ratio;
    m_slaterMatrixOld = m_slaterMatrix;
    m_slaterMatrixInverseOld = m_slaterMatrixInverse;
    m_slaterMatrixDerOld = m_slaterMatrixDer;
    m_slaterMatrixSecDerOld = m_slaterMatrixSecDer;
}

void SlaterDeterminant::resetArrays()
{
    m_positions = m_positionsOld;
    m_determinantDerivative = m_determinantDerivativeOld;
    m_determinantSecondDerivative = m_determinantSecondDerivativeOld;
    m_ratio = m_ratioOld;
    m_slaterMatrix = m_slaterMatrixOld;
    m_slaterMatrixInverse = m_slaterMatrixInverseOld;
    m_slaterMatrixDer = m_slaterMatrixDerOld;
    m_slaterMatrixSecDer = m_slaterMatrixSecDerOld;
}

void SlaterDeterminant::initializeArrays(const arma::vec positions,
                                         const arma::vec radialVector,
                                         const arma::mat distanceMatrix)
{
    m_ratio = 1;
    arma::mat positionBlock = reshape(positions, m_numberOfDimensions, m_numberOfParticles);
    m_positions = positionBlock;
    initializeSlaterMatrix();
    initializeSlaterMatrixDer();
    initializeSlaterMatrixSecDer();
    initializeSlaterMatrixInverse();
    m_determinantDerivative.zeros(m_degreesOfFreedom);
    m_determinantSecondDerivative.zeros(m_degreesOfFreedom);
    updateSlaterDeterminantDerivatives(0, m_numberOfParticles);
}

void SlaterDeterminant::updateParameters(arma::mat parameters) {}

void SlaterDeterminant::initializeSlaterMatrix()
{
    m_slaterMatrix.ones(m_numberOfParticles, m_numberOfParticlesHalf);
    for (arma::uword row = 0; row < m_numberOfParticles; row++) {
        updateSlaterMatrixRow(row);
    }
}

void SlaterDeterminant::initializeSlaterMatrixDer()
{
    m_slaterMatrixDer.zeros(m_degreesOfFreedom, m_numberOfParticlesHalf);
    for (arma::uword row = 0; row < m_degreesOfFreedom; row++) {
        updateSlaterMatrixDerRow(row);
    }
}

void SlaterDeterminant::initializeSlaterMatrixSecDer()
{
    m_slaterMatrixSecDer.zeros(m_degreesOfFreedom, m_numberOfParticlesHalf);
    for (arma::uword row = 0; row < m_degreesOfFreedom; row++) {
        updateSlaterMatrixSecDerRow(row);
    }
}

void SlaterDeterminant::initializeSlaterMatrixInverse()
{
    m_slaterMatrixInverse.zeros(m_numberOfParticlesHalf, m_numberOfParticles);
    m_slaterMatrixInverse.head_cols(m_numberOfParticlesHalf)
        = arma::inv(m_slaterMatrix.head_rows(m_numberOfParticlesHalf));
    m_slaterMatrixInverse.tail_cols(m_numberOfParticlesHalf)
        = arma::inv(m_slaterMatrix.tail_rows(m_numberOfParticlesHalf));
}

void SlaterDeterminant::updateSlaterMatrixRow(const arma::uword row)
{
    for (arma::uword col = 0; col < m_numberOfParticlesHalf; col++) {
        m_slaterMatrix(row, col) = m_basis->basisElement(col, m_positions.col(row));
    }
}

void SlaterDeterminant::updateSlaterMatrixDerRow(const arma::uword row)
{
    arma::uword particle = arma::uword(row / m_numberOfDimensions);
    arma::uword dimension = row % m_numberOfDimensions;
    for (arma::uword col = 0; col < m_numberOfParticlesHalf; col++) {
        m_slaterMatrixDer(row, col) = m_basis->basisElementDer(col,
                                                               dimension,
                                                               m_positions.col(particle));
    }
}

void SlaterDeterminant::updateSlaterMatrixSecDerRow(const arma::uword row)
{
    arma::uword particle = arma::uword(row / m_numberOfDimensions);
    arma::uword dimension = row % m_numberOfDimensions;
    for (arma::uword col = 0; col < m_numberOfParticlesHalf; col++) {
        m_slaterMatrixSecDer(row, col) = m_basis->basisElementSecDer(col,
                                                                     dimension,
                                                                     m_positions.col(particle));
    }
}

void SlaterDeterminant::updateRatio()
{
    m_ratio = arma::dot(m_slaterMatrix.row(m_particle), m_slaterMatrixInverse.col(m_particle));
}

void SlaterDeterminant::updateSlaterDeterminantDerivatives(arma::uword start, arma::uword end)
{
    for (arma::uword i = start * m_numberOfDimensions; i < end * m_numberOfDimensions; i++) {
        arma::uword particle = arma::uword(i / m_numberOfDimensions);
        m_determinantDerivative(i) = arma::dot(m_slaterMatrixDer.row(i), m_slaterMatrixInverse.col(particle));
        m_determinantSecondDerivative(i) = arma::dot(m_slaterMatrixSecDer.row(i),
                                           m_slaterMatrixInverse.col(particle));
    }
}

void SlaterDeterminant::updateSlaterMatrixInverse(arma::uword start, arma::uword end)
{
    updateRatio();
    for (arma::uword j = start; j < end; j++) {
        if (j != m_particle) {
            double S = arma::dot(m_slaterMatrix.row(m_particle), m_slaterMatrixInverse.col(j));
            m_slaterMatrixInverse.col(j) -= S * m_slaterMatrixInverse.col(m_particle) / m_ratio;
        }
    }
    m_slaterMatrixInverse.col(m_particle) /= m_ratio;
}

double SlaterDeterminant::evaluateRatio()
{
    return m_ratio * m_ratio;
}

double SlaterDeterminant::computeGradient(const arma::uword k)
{
    return m_determinantDerivative(k);
}

double SlaterDeterminant::computeLaplacian()
{
    return arma::sum(m_determinantSecondDerivative)
           - arma::dot(m_determinantDerivative, m_determinantDerivative);
}

arma::vec SlaterDeterminant::computeParameterGradient()
{
    return m_gradients;
}
