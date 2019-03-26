#include "slaterdeterminant.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../Basis/basis.h"

SlaterDeterminant::SlaterDeterminant(System* system) :
        WaveFunction(system) {
    m_numberOfFreeDimensions            = m_system->getNumberOfFreeDimensions();
    m_freeDimensionsHalf                = m_numberOfFreeDimensions/2;
    m_numberOfParticles                 = m_system->getNumberOfParticles();
    //m_numberOfOrbitals                  = m_system->getBasis()->getNumberOfOrbitals();
    m_numberOfDimensions                = m_system->getNumberOfDimensions();
    m_numberOfParticlesHalf             = m_numberOfParticles/2;
    m_maxNumberOfParametersPerElement   = m_system->getMaxNumberOfParametersPerElement();
}



void SlaterDeterminant::updateArrays(const Eigen::VectorXd positions, const int pRand) {
    //Update old arrays
    m_positionsOld              = m_positions;
    m_determinantDerivativeOld  = m_determinantDerivative;

    //Update new arrays
    m_positions = positions;
    int particle = int(pRand/m_numberOfDimensions);

    m_probabilityRatioOld = m_probabilityRatio;
    m_slaterMatrixOld     = m_slaterMatrix;
    m_slaterMatrixInverseOld = m_slaterMatrixInverse;

    if(particle < m_numberOfParticlesHalf) {
        //m_slaterMatrixUpOld                 = m_slaterMatrixUp;
        m_slaterMatrixUpDerOld              = m_slaterMatrixUpDer;
        m_slaterMatrixUpSecDerOld           = m_slaterMatrixUpSecDer;
        //m_slaterMatrixUpInverseOld          = m_slaterMatrixUpInverse;

        m_slaterMatrix.row(particle)    = updateSlaterMatrixRow(m_positions.head(m_freeDimensionsHalf), particle);

        for(int i=0; i<m_numberOfDimensions; i++) {
            int k = m_numberOfDimensions * particle + i;
            m_slaterMatrixUpDer.row(k)   = updateSlaterMatrixDerRow(m_positions.head(m_freeDimensionsHalf), k);
        }

        double R = 0;
        for(int j=0; j<m_numberOfParticlesHalf; j++) {
            R += m_slaterMatrix(particle, j) * m_slaterMatrixInverse(j, particle);
        }
        for(int j=0; j<m_numberOfParticlesHalf; j++) {
            double S = 0;
            for(int l=0; l<m_numberOfParticlesHalf; l++) {
                S += m_slaterMatrix(particle, l) * m_slaterMatrixInverse(l,j);
            }
            if(j != particle) {
                m_slaterMatrixInverse.col(j) -= S * m_slaterMatrixInverse.col(particle) / R;
            }
        }
        m_slaterMatrixInverse.col(particle) /= R;
        for(int i=0; i<m_freeDimensionsHalf; i++) {
            m_determinantDerivative(i) = 0;
            m_determinantSecondDerivative(i) = 0;
            int l = int(i/m_numberOfDimensions);
            for(int j=0; j<m_numberOfParticlesHalf; j++) {
                m_determinantDerivative(i)       += m_slaterMatrixUpDer(i,j)    * m_slaterMatrixInverse(j,l);
                m_determinantSecondDerivative(i) += m_slaterMatrixUpSecDer(i,j) * m_slaterMatrixInverse(j,l);
            }
        }
        m_probabilityRatio = R*R;
    }
    else {
        //m_slaterMatrixDnOld                 = m_slaterMatrixDn;
        m_slaterMatrixDnDerOld              = m_slaterMatrixDnDer;
        m_slaterMatrixDnSecDerOld           = m_slaterMatrixDnSecDer;
        //m_slaterMatrixDnInverseOld          = m_slaterMatrixDnInverse;

        int particle2 = particle - m_numberOfParticlesHalf;
        m_slaterMatrix.row(particle)    = updateSlaterMatrixRow(m_positions.tail(m_freeDimensionsHalf), particle2);
        for(int i=0; i<m_numberOfDimensions; i++) {
            int k = m_numberOfDimensions * particle2 + i;
            m_slaterMatrixDnDer.row(k)   = updateSlaterMatrixDerRow(m_positions.tail(m_freeDimensionsHalf), k);
        }
        double R = 0;
        for(int j=0; j<m_numberOfParticlesHalf; j++) {
            R += m_slaterMatrix(particle, j) * m_slaterMatrixInverse(j, particle);
        }
        for(int j=m_numberOfParticlesHalf; j<m_numberOfParticles; j++) {
            double S = 0;
            for(int l=0; l<m_numberOfParticlesHalf; l++) {
                S += m_slaterMatrix(particle, l) * m_slaterMatrixInverse(l,j);
            }
            if(j != particle) {
                m_slaterMatrixInverse.col(j) -= S * m_slaterMatrixInverse.col(particle) / R;
            }
        }
        m_slaterMatrixInverse.col(particle) /= R;
        for(int i=m_freeDimensionsHalf; i<m_numberOfFreeDimensions; i++) {
            m_determinantDerivative(i) = 0;
            m_determinantSecondDerivative(i) = 0;
            int k = i - m_freeDimensionsHalf;
            int l = int(k/m_numberOfDimensions) + m_numberOfParticlesHalf;
            for(int j=0; j<m_numberOfParticlesHalf; j++) {
                m_determinantDerivative(i)       += m_slaterMatrixDnDer(k,j)    * m_slaterMatrixInverse(j,l);
                m_determinantSecondDerivative(i) += m_slaterMatrixDnSecDer(k,j) * m_slaterMatrixInverse(j,l);
            }
        }
        m_probabilityRatio = R*R;
    }
    //std::cout << "gjgj" << std::endl;
}

void SlaterDeterminant::resetArrays(int pRand) {
    m_positions             = m_positionsOld;
    m_determinantDerivative = m_determinantDerivativeOld;
    m_determinantSecondDerivative = m_determinantSecondDerivativeOld;
    m_probabilityRatio      = m_probabilityRatioOld;
    m_slaterMatrix          = m_slaterMatrixOld;
    m_slaterMatrixInverse   = m_slaterMatrixInverseOld;

    if(pRand < m_freeDimensionsHalf) {
        //m_slaterMatrixUp            = m_slaterMatrixUpOld;
        m_slaterMatrixUpDer         = m_slaterMatrixUpDerOld;
        //m_slaterMatrixUpInverse     = m_slaterMatrixUpInverseOld;
    }
    else {
        //m_slaterMatrixDn            = m_slaterMatrixDnOld;
        m_slaterMatrixDnDer         = m_slaterMatrixDnDerOld;
        //m_slaterMatrixDnInverse     = m_slaterMatrixDnInverseOld;
    }
}

void SlaterDeterminant::initializeArrays(const Eigen::VectorXd positions) {

    m_positions                 = positions;
    m_listOfStates              = m_system->getBasis()->generateListOfStates();

    //m_slaterMatrixUp            = initializeSlaterMatrix(m_positions.head(m_numberOfFreeDimensions/2));
    //m_slaterMatrixDn            = initializeSlaterMatrix(m_positions.tail(m_numberOfFreeDimensions/2));
    initializeSlaterMatrix();

    //m_slaterMatrixUpOld         = m_slaterMatrixUp;
    //m_slaterMatrixDnOld         = m_slaterMatrixDn;
    m_slaterMatrixOld           = m_slaterMatrix;

    //m_slaterMatrixUpInverse     = m_slaterMatrix.topRows(m_numberOfParticlesHalf).inverse();
    //m_slaterMatrixDnInverse     = m_slaterMatrix.bottomRows(m_numberOfParticlesHalf).inverse();
    m_slaterMatrixInverse = Eigen::MatrixXd::Zero(m_numberOfParticlesHalf, m_numberOfParticles);
    m_slaterMatrixInverse.leftCols (m_numberOfParticlesHalf) = m_slaterMatrix.topRows   (m_numberOfParticlesHalf).inverse();
    m_slaterMatrixInverse.rightCols(m_numberOfParticlesHalf) = m_slaterMatrix.bottomRows(m_numberOfParticlesHalf).inverse();

    //m_slaterMatrixUpInverseOld  = m_slaterMatrixUpInverse;
    //m_slaterMatrixDnInverseOld  = m_slaterMatrixDnInverse;
    m_slaterMatrixInverseOld = m_slaterMatrixInverse;

    m_slaterMatrixUpDer         = initializeSlaterMatrixDer(m_positions.head(m_numberOfFreeDimensions/2));
    m_slaterMatrixDnDer         = initializeSlaterMatrixDer(m_positions.tail(m_numberOfFreeDimensions/2));
    m_slaterMatrixUpDerOld      = m_slaterMatrixUpDer;
    m_slaterMatrixDnDerOld      = m_slaterMatrixDnDer;

    m_slaterMatrixUpSecDer      = initializeSlaterMatrixSecDer(m_positions.head(m_numberOfFreeDimensions/2));
    m_slaterMatrixDnSecDer      = initializeSlaterMatrixSecDer(m_positions.tail(m_numberOfFreeDimensions/2));
    m_slaterMatrixUpSecDerOld   = m_slaterMatrixUpSecDer;
    m_slaterMatrixDnSecDerOld   = m_slaterMatrixDnSecDer;

    m_determinantDerivative     = Eigen::VectorXd::Zero(m_numberOfFreeDimensions);
    m_determinantSecondDerivative  = Eigen::VectorXd::Zero(m_numberOfFreeDimensions);
    for(int i=0; i<m_freeDimensionsHalf; i++) {
        int k = i + m_freeDimensionsHalf;
        for(int j=0; j<m_numberOfParticlesHalf; j++) {
            m_determinantDerivative(i) += m_slaterMatrixUpDer(i,j) * m_slaterMatrixInverse(j,int(i/m_numberOfDimensions));
            m_determinantDerivative(k) += m_slaterMatrixDnDer(i,j) * m_slaterMatrixInverse(j,int(i/m_numberOfDimensions)+m_numberOfParticlesHalf);
            m_determinantSecondDerivative(i) += m_slaterMatrixUpSecDer(i,j) * m_slaterMatrixInverse(j,int(i/m_numberOfDimensions));
            m_determinantSecondDerivative(k) += m_slaterMatrixDnSecDer(i,j) * m_slaterMatrixInverse(j,int(i/m_numberOfDimensions)+m_numberOfParticlesHalf);
        }
    }
    m_determinantDerivativeOld       = m_determinantDerivative;
    m_determinantSecondDerivativeOld = m_determinantSecondDerivative;
    m_probabilityRatio               = 1;
    m_probabilityRatioOld            = m_probabilityRatio;
}

void SlaterDeterminant::updateParameters(const Eigen::MatrixXd parameters, const int elementNumber) {
    m_elementNumber = elementNumber;
}

Eigen::VectorXd SlaterDeterminant::updateSlaterMatrixDerRow(const Eigen::VectorXd positions, const int k) {
    //Update row of dA
    Eigen::VectorXd dA = Eigen::VectorXd::Zero(m_numberOfParticlesHalf);
    int particle  = int(k/m_numberOfDimensions);
    int dimension = k%m_numberOfDimensions;
    // Find matrix
    for(int i=0; i<m_numberOfParticlesHalf; i++) {
        dA(i) = m_system->getBasis()->evaluateDerivative(positions(k), int(m_listOfStates(i, dimension)));
        for(int j=0; j<m_numberOfDimensions; j++) {
            int m = m_numberOfDimensions * particle + j;
            if(m != k) {
                dA(i) *= m_system->getBasis()->evaluate(positions(m), int(m_listOfStates(i, j)));
            }
        }
    }
    return dA.transpose();
}

Eigen::MatrixXd SlaterDeterminant::initializeSlaterMatrixDer(const Eigen::VectorXd positions) {
    //Initialize the entire dA matrix
    Eigen::MatrixXd dA = Eigen::MatrixXd::Zero(m_freeDimensionsHalf, m_numberOfParticlesHalf);
    for(int k=0; k<m_freeDimensionsHalf; k++) {
        dA.row(k) = updateSlaterMatrixDerRow(positions, k);
    }
    return dA;
}

Eigen::VectorXd SlaterDeterminant::updateSlaterMatrixSecDerRow(const Eigen::VectorXd positions, const int k) {
    //Update row of dA
    Eigen::VectorXd d2A = Eigen::VectorXd::Zero(m_numberOfParticlesHalf);
    int particle  = int(k/m_numberOfDimensions);
    int dimension = k%m_numberOfDimensions;
    // Find matrix
    for(int i=0; i<m_numberOfParticlesHalf; i++) {
        d2A(i) = m_system->getBasis()->evaluateSecondDerivative(positions(k), int(m_listOfStates(i, dimension)));
        for(int j=0; j<m_numberOfDimensions; j++) {
            int m = m_numberOfDimensions * particle + j;
            if(m != k) {
                d2A(i) *= m_system->getBasis()->evaluate(positions(m), int(m_listOfStates(i, j)));
            }
        }
    }
    return d2A;
}

Eigen::MatrixXd SlaterDeterminant::initializeSlaterMatrixSecDer(const Eigen::VectorXd positions) {
    //Initialize the entire dA matrix
    Eigen::MatrixXd d2A = Eigen::MatrixXd::Zero(m_freeDimensionsHalf, m_numberOfParticlesHalf);
    for(int k=0; k<m_freeDimensionsHalf; k++) {
        d2A.row(k) = updateSlaterMatrixSecDerRow(positions, k);
    }
    return d2A;
}

double SlaterDeterminant::updateSlaterMatrixElement(const Eigen::VectorXd positions, const int i, const int j) {
    // Updates an element in A-matrix
    double element = 1;
    for(int k=0; k<m_numberOfDimensions; k++) {
        element *= m_system->getBasis()->evaluate(positions(m_numberOfDimensions*i+k), int(m_listOfStates(j,k)));
    }
    return element;
}

Eigen::VectorXd SlaterDeterminant::updateSlaterMatrixRow(const Eigen::VectorXd positions, const int i) {
    // Updates a row in A-matrix
    Eigen::VectorXd A = Eigen::VectorXd::Ones(m_numberOfParticlesHalf);
    for(int j=0; j<m_numberOfParticlesHalf; j++) {
        A(j) = updateSlaterMatrixElement(positions, i, j);
    }
    return A;
}

void SlaterDeterminant::initializeSlaterMatrix() {
    // Update the entire matrix
    m_slaterMatrix = Eigen::MatrixXd::Ones(m_numberOfParticles, m_numberOfParticlesHalf);
    for(int j=0; j<m_numberOfParticles; j++) {
        m_slaterMatrix.row(j) = updateSlaterMatrixRow(m_positions, j);
    }
}

double SlaterDeterminant::evaluateRatio() {
    return m_probabilityRatio;
}

double SlaterDeterminant::computeFirstDerivative(const int k) {
    return m_determinantDerivative(k);
}

double SlaterDeterminant::computeSecondDerivative() {
    return m_determinantSecondDerivative.cwiseAbs2().sum() - double(m_determinantDerivative.transpose() * m_determinantDerivative);
}

Eigen::VectorXd SlaterDeterminant::computeFirstEnergyDerivative(const int k) {
    return Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);
}

Eigen::VectorXd SlaterDeterminant::computeSecondEnergyDerivative() {
    return Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);
}
