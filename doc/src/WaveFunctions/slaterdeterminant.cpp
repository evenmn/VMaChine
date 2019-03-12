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
    m_numberOfOrbitals                  = m_system->getBasis()->getNumberOfOrbitals();
    m_numberOfDimensions                = m_system->getNumberOfDimensions();
    m_numberOfParticlesHalf             = m_numberOfParticles/2;
    m_maxNumberOfParametersPerElement   = m_system->getMaxNumberOfParametersPerElement();
    m_omega                             = m_system->getFrequency();
}

void SlaterDeterminant::updateArrays(const Eigen::VectorXd positions, const int pRand) {
    //Update old arrays
    m_oldPositions          = m_positions;
    m_diffOld               = m_diff;

    //Update new arrays
    m_positions = positions;
    int particle = int(pRand/m_numberOfDimensions);

    m_oldRatio = m_ratio;

    if(particle < m_numberOfParticlesHalf) {
        m_D_upOld               = m_D_up;
        m_dD_upOld              = m_dD_up;
        m_D_up_old_inv          = m_D_up_inv;
        m_old_det_up            = m_det_up;

        m_D_up.row(particle)    = updateRow(m_positions.head(m_freeDimensionsHalf), particle);
        for(int i=0; i<m_numberOfDimensions; i++) {
            int k = m_numberOfDimensions * particle + i;
            m_dD_up.row(k)   = dA_row(m_positions.head(m_freeDimensionsHalf), k);
        }
        double R = 0;
        for(int j=0; j<m_numberOfParticlesHalf; j++) {
            R += m_D_up(particle, j) * m_D_up_old_inv(j, particle);
        }
        for(int j=0; j<m_numberOfParticlesHalf; j++) {
            double S = 0;
            for(int l=0; l<m_numberOfParticlesHalf; l++) {
                S += m_D_up(particle, l) * m_D_up_inv(l,j);
            }
            if(j != particle) {
                m_D_up_inv.col(j) -= S * m_D_up_inv.col(particle) / R;
            }
        }
        m_D_up_inv.col(particle) /= R;
        for(int i=0; i<m_freeDimensionsHalf; i++) {
            m_diff(i) = 0;
            for(int j=0; j<m_numberOfParticlesHalf; j++) {
                m_diff(i) += m_dD_up(i,j) * m_D_up_inv(j,int(i/2));
            }
        }

        m_det_up = m_D_up.determinant();
        m_ratio  = m_det_up * m_det_up / (m_old_det_up * m_old_det_up);
    }
    else {
        m_D_dnOld               = m_D_dn;
        m_dD_dnOld              = m_dD_dn;
        m_D_dn_old_inv          = m_D_dn_inv;
        m_old_det_dn            = m_det_dn;

        int particle2 = particle - m_numberOfParticlesHalf;
        m_D_dn.row(particle2)   = updateRow(m_positions.tail(m_freeDimensionsHalf), particle2);
        for(int i=0; i<m_numberOfDimensions; i++) {
            int k = m_numberOfDimensions * particle2 + i;
            m_dD_dn.row(k)   = dA_row(m_positions.tail(m_freeDimensionsHalf), k);
        }
        double R = 0;
        for(int j=0; j<m_numberOfParticlesHalf; j++) {
            R += m_D_dn(particle2, j) * m_D_dn_old_inv(j, particle2);
        }
        for(int j=0; j<m_numberOfParticlesHalf; j++) {
            double S = 0;
            for(int l=0; l<m_numberOfParticlesHalf; l++) {
                S += m_D_dn(particle2, l) * m_D_dn_inv(l,j);
            }
            if(j != particle2) {
                m_D_dn_inv.col(j) -= S * m_D_dn_inv.col(particle2) / R;
            }
        }
        m_D_dn_inv.col(particle2) /= R;
        for(int i=m_freeDimensionsHalf; i<m_numberOfFreeDimensions; i++) {
            m_diff(i) = 0;
            int k = i-m_freeDimensionsHalf;
            for(int j=0; j<m_numberOfParticlesHalf; j++) {
                m_diff(i) += m_dD_dn(k,j) * m_D_dn_inv(j,int(k/2));
            }
        }
        m_det_dn = m_D_dn.determinant();
        m_ratio  = m_det_dn * m_det_dn / (m_old_det_dn * m_old_det_dn);
    }
}

void SlaterDeterminant::resetArrays() {
    m_positions = m_oldPositions;
    m_D_up      = m_D_upOld;
    m_D_dn      = m_D_dnOld;
    m_D_up_inv  = m_D_up_old_inv;
    m_D_dn_inv  = m_D_dn_old_inv;
    m_dD_up     = m_dD_upOld;
    m_dD_dn     = m_dD_dnOld;
    m_diff      = m_diffOld;
    m_det_up    = m_old_det_up;
    m_det_dn    = m_old_det_dn;
    m_ratio     = m_oldRatio;
}

void SlaterDeterminant::initializeArrays(const Eigen::VectorXd positions) {
    m_positions     = positions;
    m_order         = list();

    m_D_up          = updateMatrix(m_positions.head(m_numberOfFreeDimensions/2));
    m_D_dn          = updateMatrix(m_positions.tail(m_numberOfFreeDimensions/2));
    m_D_upOld       = m_D_up;
    m_D_dnOld       = m_D_dn;

    m_D_up_inv      = m_D_up.inverse();
    m_D_dn_inv      = m_D_dn.inverse();
    m_D_up_old_inv  = m_D_up_inv;
    m_D_dn_old_inv  = m_D_dn_inv;

    m_dD_up         = dA_matrix(m_positions.head(m_numberOfFreeDimensions/2));
    m_dD_dn         = dA_matrix(m_positions.tail(m_numberOfFreeDimensions/2));
    m_dD_upOld      = m_dD_up;
    m_dD_dnOld      = m_dD_dn;

    m_det_up        = m_D_up.determinant();
    m_det_dn        = m_D_dn.determinant();
    m_old_det_up    = m_det_up;
    m_old_det_dn    = m_det_dn;

    m_diff  = Eigen::VectorXd::Zero(m_numberOfFreeDimensions);
    for(int i=0; i<m_numberOfFreeDimensions/2; i++) {
        for(int j=0; j<m_numberOfParticlesHalf; j++) {
            m_diff(i) += m_dD_up(i,j) * m_D_up_inv(j,int(i/2));
            m_diff(i+m_numberOfFreeDimensions/2) += m_dD_dn(i,j) * m_D_dn_inv(j,int(i/2));
        }
    }
    m_diffOld = m_diff;

    m_ratio = 1;
}

void SlaterDeterminant::updateParameters(const Eigen::MatrixXd parameters, const int elementNumber) {
    m_elementNumber = elementNumber;
}

Eigen::MatrixXd SlaterDeterminant::list() {
    // Returns the index list used in Slater
    // For instance (0,0), (1,0), (0,1) for 6P in 2D
    //              (0,0,0), (1,0,0), (0,1,0), (0,0,1) for 6P in 3D etc..
    Eigen::MatrixXd order = Eigen::MatrixXd::Zero(m_numberOfParticlesHalf, m_numberOfDimensions);
    int counter = 0;
    // Two dimensions
    if (m_numberOfDimensions == 2) {
        for(int i=0; i<m_numberOfOrbitals; i++) {
            for(int s=i; s<m_numberOfOrbitals; s++) {
                int j = s - i;
                order(counter,1) = i;
                order(counter,0) = j;
                counter += 1;
            }
        }
    }
    // Three dimensions
    else if (m_numberOfDimensions == 3) {
        for(int i=0; i<m_numberOfOrbitals; i++) {
            for(int j=0; j<m_numberOfOrbitals; j++) {
                for(int s=i+j; s<m_numberOfOrbitals; s++) {
                    int k = s - i - j;
                    order(counter,0) = i;
                    order(counter,1) = j;
                    order(counter,2) = k;
                    counter += 1;
                }
            }
        }
    }

    else {
        std::cout << "Number of dimensions should either be 2 or 3" << std::endl;
        exit(0);
    }
    return order;
}

Eigen::VectorXd SlaterDeterminant::dA_row(const Eigen::VectorXd positions, const int k) {
    //Update row of dA
    Eigen::VectorXd dA = Eigen::VectorXd::Zero(m_numberOfParticlesHalf);
    int particle  = int(k/m_numberOfDimensions);
    int dimension = k%m_numberOfDimensions;
    // Find matrix
    for(int i=0; i<m_numberOfParticlesHalf; i++) {
        dA(i) = m_system->getBasis()->evaluateDerivative(sqrt(m_omega) * positions(k), int(m_order(i, dimension)));
        for(int j=0; j<m_numberOfDimensions; j++) {
            int m = m_numberOfDimensions * particle + j;
            if(m != k) {
                dA(i) *= m_system->getBasis()->evaluate(sqrt(m_omega) * positions(m), int(m_order(i, j)));
            }
        }
    }
    return dA.transpose();
}

Eigen::MatrixXd SlaterDeterminant::dA_matrix(const Eigen::VectorXd positions) {
    //Initialize the entire dA matrix
    Eigen::MatrixXd dA = Eigen::MatrixXd::Zero(m_numberOfParticles, m_numberOfParticlesHalf);
    for(int k=0; k<m_numberOfParticlesHalf; k++) {
        dA.row(k) = dA_row(positions, k);
    }
    return dA;
}

double SlaterDeterminant::updateElement(const Eigen::VectorXd positions, const int i, const int j) {
    // Updates an element in A-matrix
    double element = 1;
    for(int k=0; k<m_numberOfDimensions; k++) {
        element *= m_system->getBasis()->evaluate(sqrt(m_omega) * positions(m_numberOfDimensions*i+k), int(m_order(j,k)));
    }
    return element;
}

Eigen::VectorXd SlaterDeterminant::updateRow(const Eigen::VectorXd positions, const int i) {
    // Updates a row in A-matrix
    Eigen::VectorXd A = Eigen::VectorXd::Ones(m_numberOfParticlesHalf);
    for(int j=0; j<m_numberOfParticlesHalf; j++) {
        A(j) = updateElement(positions, i, j);
    }
    return A;
}

Eigen::MatrixXd SlaterDeterminant::updateMatrix(const Eigen::VectorXd positions) {
    // Update the entire matrix
    Eigen::MatrixXd A = Eigen::MatrixXd::Ones(m_numberOfParticlesHalf, m_numberOfParticlesHalf);
    for(int j=0; j<m_numberOfParticlesHalf; j++) {
        A.row(j) = updateRow(positions, j);
    }
    return A;
}

double SlaterDeterminant::evaluateRatio() {
    return m_ratio;
}

double SlaterDeterminant::computeFirstDerivative(const int k) {
    return m_diff(k);
}

double SlaterDeterminant::computeSecondDerivative() {
    return -double(m_diff.transpose() * m_diff);
}

Eigen::VectorXd SlaterDeterminant::computeFirstEnergyDerivative(const int k) {
    return Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);
}

Eigen::VectorXd SlaterDeterminant::computeSecondEnergyDerivative() {
    return Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);
}
