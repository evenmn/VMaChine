#include "hartreefock.h"
#include "hermite.h"
#include "hermiteexpansion.h"
#include "hermitespin.h"
#include "hydrogenorbital.h"
#include "none.h"
#include "../system.h"
#include <iostream>
#include <fstream>

HartreeFock::HartreeFock(System *system, Basis *basis)  :
    Basis(system) {
    m_system                = system;
    m_numberOfParticles     = m_system->getNumberOfParticles();
    m_numberOfDimensions    = m_system->getNumberOfDimensions();
    m_omega                 = m_system->getFrequency();
    m_omegaSqrt             = sqrt(m_omega);
    m_path                  = m_system->getPath();
    m_basis                 = basis;
    readCoefficientFile();
    Basis::numberOfOrbitals();
    Basis::generateListOfStates();
}

std::string HartreeFock::generateFileName() {
    std::string fileName = m_path;
    fileName += "int1/";
    fileName += "hartree-fock/";
    fileName += "harmonicoscillator/";
    fileName += std::to_string(m_numberOfDimensions) + "D/";
    fileName += std::to_string(m_numberOfParticles) + "P/";
    fileName += std::to_string(m_omega) + "w/";
    fileName += "coeffs.dat";
    return fileName;
}

void HartreeFock::readCoefficientFile() {
    std::string fileName = generateFileName();
    std::ifstream inFile(generateFileName().c_str(), std::ios::in);
    m_basisSize     = Basis::fileLength(fileName);
    m_coefficients  = Eigen::MatrixXd::Zero(int(m_numberOfParticles/2), m_basisSize);
    if (!inFile.is_open()) {
        std::cout << "file not found";
        MPI_Finalize();
        exit(0);
    }
    else {
        double value;
        int counter = 0;
        while (inFile >> value) {
            m_coefficients(int(counter/m_basisSize), counter % m_basisSize) = value;
            counter += 1;
        }
    }
}
/*
void HartreeFock::numberOfOrbitals() {
    //Number of closed-shell orbitals
    int counter = 0;
    while(true) {
        int orb = Basis::binomial(counter, m_numberOfDimensions);
        if(orb == m_basisSize) {
            m_numberOfOrbitals = counter+1;
            break;
        }
        else if(orb > m_basisSize) {
            std::cout << "Basis size must correspond to a closed shell..." << std::endl;
            MPI_Finalize();
            exit(0);
        }
        counter += 1;
    }
}

void HartreeFock::generateListOfStates(int orbitals) {
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
        exit(0);
    }
}
*/

double HartreeFock::evaluate(double x, int n) {
    //Hermite polynomial of n'th degree
    double sum = 0;
    for(int lambda=0; lambda<m_basisSize; lambda++) {
        sum += m_coefficients(n, lambda) * m_basis->evaluate(x, lambda);
    }
    return sum;
}

double HartreeFock::evaluateDerivative(double x, int n) {
    //First derivative of Hermite polynomial of n'th degree
    double sum = 0;
    for(int lambda=0; lambda<m_basisSize; lambda++) {
        sum += m_coefficients(n, lambda) * m_basis->evaluateDerivative(x, lambda);
    }
    return sum;
}

double HartreeFock::evaluateSecondDerivative(const double x, const int n) {
    //Second derivative of Hermite polynomial of n'th degree
    double sum = 0;
    for(int lambda=0; lambda<m_basisSize; lambda++) {
        sum += m_coefficients(n, lambda) * m_basis->evaluateSecondDerivative(x, lambda);
    }
    return sum;
}

double HartreeFock::basisElement(const int n, Eigen::VectorXd positions) {
    double prod = 1;
    for(int i=0; i<m_numberOfDimensions; i++) {
        prod *= evaluate(positions(i), int(m_listOfStates(n, i)));
    }
    return prod;
}

double HartreeFock::basisElementDer(const int n, const int i, Eigen::VectorXd positions) {
    // i is the dimension we are derivating with respect to
    double prod = evaluateDerivative(positions(i), m_listOfStates(n, i));
    for(int j=0; j<m_numberOfDimensions; j++) {
        if(i != j) {
            prod *= evaluate(positions(j), m_listOfStates(n, j));
        }
    }
    return prod;
}

double HartreeFock::basisElementSecDer(const int n, const int i, Eigen::VectorXd positions) {
    // i is the dimension we are derivating with respect to
    double prod = evaluateSecondDerivative(positions(i), m_listOfStates(n, i));
    for(int j=0; j<m_numberOfDimensions; j++) {
        if(i != j) {
            prod *= evaluate(positions(j), m_listOfStates(n, j));
        }
    }
    return prod;
}
