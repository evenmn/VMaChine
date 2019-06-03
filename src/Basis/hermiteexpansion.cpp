#include "hermiteexpansion.h"
#include "hermite.h"
#include "../system.h"
#include <iostream>
#include <fstream>

HermiteExpansion::HermiteExpansion(System *system)  :
    Basis(system) {
    m_system                = system;
    m_numberOfParticles     = m_system->getNumberOfParticles();
    m_numberOfDimensions    = m_system->getNumberOfDimensions();
    m_omega                 = m_system->getFrequency();
    m_path                  = m_system->getPath();
    m_basis                 = new Hermite(system);
    readCoefficientFile();
    numberOfOrbitals();
    generateListOfStates(m_numberOfOrbitals);
}

std::string HermiteExpansion::generateFileName() {
    std::string fileName = m_path;
    fileName += "int1/";
    fileName += "doublewell/";
    fileName += std::to_string(m_omega) + "w/";
    fileName += "coeffs.dat";
    return fileName;
}

void HermiteExpansion::readCoefficientFile() {
    std::string fileName = generateFileName();
    m_basisSize     = Basis::fileLength(fileName);
    m_coefficients  = Eigen::MatrixXd::Zero(m_basisSize, m_basisSize);
    Basis::writeFileContentIntoEigenMatrix(fileName, m_coefficients);
}

void HermiteExpansion::numberOfOrbitals() {
    //Number of closed-shell orbitals
    int counter = 0;
    while(true) {
        int orb = 2 * Basis::binomial(counter, m_numberOfDimensions);
        if(orb == m_numberOfParticles) {
            m_numberOfOrbitals = counter+1;
            break;
        }
        else if(orb > m_numberOfParticles) {
            std::cout << "This program supports closed-shells only. Please choose a "
                         "number of particles such that the orbital is full" << std::endl;
            MPI_Finalize();
            exit(0);
        }
        counter += 1;
    }
}

void HermiteExpansion::generateListOfStates(int orbitals) {
    // Returns the index list used in Slater
    // For instance (0,0), (1,0), (0,1) for 6P in 2D
    //              (0,0,0), (1,0,0), (0,1,0), (0,0,1) for 8P in 3D etc..
    int numberOfStates = Basis::binomial(orbitals-1, m_numberOfDimensions);
    m_listOfStates = Eigen::MatrixXi::Zero(numberOfStates, m_numberOfDimensions);
    int counter = 0;
    // Two dimensions
    if (m_numberOfDimensions == 2) {
        for(int i=0; i<orbitals; i++) {
            for(int j=0; j<i+1; j++) {
                m_listOfStates(counter,0) = i-j;
                m_listOfStates(counter,1) = j;
                counter += 1;
            }
        }
    }
    // Three dimensions
    else if (m_numberOfDimensions == 3) {
        for(int i=0; i<orbitals; i++) {
            for(int j=0; j<i+1; j++) {
                for(int k=0; k<i-j+1; k++) {
                    m_listOfStates(counter,0) = i-j-k;
                    m_listOfStates(counter,1) = j;
                    m_listOfStates(counter,2) = k;
                    counter += 1;
                }
            }
        }
    }
    else {
        std::cout << "Number of dimensions should be either 2 or 3" << std::endl;
        MPI_Finalize();
        exit(0);
    }
}

double HermiteExpansion::evaluate(double x, int n) {
    //Hermite polynomial of n'th degree
    return m_basis->evaluate(x, n);
}

double HermiteExpansion::evaluateDerivative(double x, int n) {
    //First derivative of Hermite polynomial of n'th degree
    return m_basis->evaluateDerivative(x, n);
}

double HermiteExpansion::evaluateSecondDerivative(const double x, const int n) {
    //Second derivative of Hermite polynomial of n'th degree
    return m_basis->evaluateSecondDerivative(x, n);
}

double HermiteExpansion::evaluateXDir(double positionx, int nx) {
    double sum = 0;
    for(int lambda=0; lambda<m_basisSize; lambda++) {
        sum += m_coefficients(lambda, nx) * evaluate(positionx, lambda);
    }
    return sum;
}

double HermiteExpansion::evaluateXDirDerivative(double positionx, int nx) {
    double sum = 0;
    for(int lambda=0; lambda<m_basisSize; lambda++) {
        sum += m_coefficients(lambda, nx) * evaluateDerivative(positionx, lambda);
    }
    return sum;
}

double HermiteExpansion::evaluateXDirSecondDerivative(double positionx, int nx) {
    double sum = 0;
    for(int lambda=0; lambda<m_basisSize; lambda++) {
        sum += m_coefficients(lambda, nx) * evaluateSecondDerivative(positionx, lambda);
    }
    return sum;
}

double HermiteExpansion::basisElement(const int n, Eigen::VectorXd positions) {
    double prod = evaluateXDir(positions(0), m_listOfStates(n, 0));
    for(int i=1; i<m_numberOfDimensions; i++) {
        prod *= evaluate(positions(i), m_listOfStates(n, i));
    }
    return prod;
}

double HermiteExpansion::basisElementDer(const int n, const int i, Eigen::VectorXd positions) {
    // i is the dimension we are derivating with respect to
    double prod;
    if(i == 0) {
        prod = evaluateXDirDerivative(positions(0), m_listOfStates(n, 0));
        for(int j=1; j<m_numberOfDimensions; j++) {
            prod *= evaluate(positions(j), m_listOfStates(n, j));
        }
    }
    else {
        prod = evaluateXDir(positions(0), m_listOfStates(n, 0));
        prod *= evaluateDerivative(positions(i), m_listOfStates(n, i));
        for(int j=1; j<m_numberOfDimensions; j++) {
            if(j != i) {
                prod *= evaluate(positions(j), m_listOfStates(n, j));
            }
        }
    }
    return prod;
}

double HermiteExpansion::basisElementSecDer(const int n, const int i, Eigen::VectorXd positions) {
    // i is the dimension we are derivating with respect to
    double prod;
    if(i == 0) {
        prod = evaluateXDirSecondDerivative(positions(0), m_listOfStates(n, 0));
        for(int j=1; j<m_numberOfDimensions; j++) {
            prod *= evaluate(positions(j), m_listOfStates(n, j));
        }
    }
    else {
        prod = evaluateXDir(positions(0), m_listOfStates(n, 0));
        prod *= evaluateSecondDerivative(positions(i), m_listOfStates(n, i));
        for(int j=1; j<m_numberOfDimensions; j++) {
            if(j != i) {
                prod *= evaluate(positions(j), m_listOfStates(n, j));
            }
        }
    }
    return prod;
}