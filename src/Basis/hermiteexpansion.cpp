#include "hermiteexpansion.h"
#include "hermite.h"
#include "../system.h"
#include "../Hamiltonians/hamiltonian.h"
#include <iostream>
#include <fstream>

HermiteExpansion::HermiteExpansion(System *system)  :
    Basis(system) {
    m_system                = system;
    m_numberOfParticles     = m_system->getNumberOfParticles();
    m_numberOfDimensions    = m_system->getNumberOfDimensions();
    m_numberOfSources       = m_system->getHamiltonian()->getNumberOfSources();
    m_omega                 = m_system->getFrequency();
    m_path                  = m_system->getPath();
    m_basis                 = new Hermite(system);
    readCoefficientFile();
    m_listOfStates = Basis::generateListOfStates(m_numberOfSources);

    /*
    int Ngrid = 1000;
    Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(Ngrid,-10,10);
    m_dim = 0;
    Eigen::VectorXd y = Eigen::VectorXd::Zero(Ngrid);
    for(int i=0; i<Ngrid; i++) {
        y(i) = evaluate(x(i), 1) * exp(-0.5*x(i)*x(i));
    }
    std::ofstream file;
    file.open(m_path + "test.dat");
    file << y << std::endl;
    */
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

double HermiteExpansion::evaluate(double x, int n) {
    //Hermite polynomial of n'th degree
    if(m_dim == 0) {
        double sum = 0;
        for(int lambda=0; lambda<m_basisSize; lambda++) {
            sum += m_coefficients(lambda, n) * m_basis->evaluate(x, lambda);
        }
        return sum;
    }
    else {
        return m_basis->evaluate(x, n);
    }
}

double HermiteExpansion::evaluateDerivative(double x, int n) {
    //First derivative of Hermite polynomial of n'th degree
    if(m_dim == 0) {
        double sum = 0;
        for(int lambda=0; lambda<m_basisSize; lambda++) {
            sum += m_coefficients(lambda, n) * m_basis->evaluateDerivative(x, lambda);
        }
        return sum;
    }
    else {
        return m_basis->evaluateDerivative(x, n);
    }
}

double HermiteExpansion::evaluateSecondDerivative(const double x, const int n) {
    //Second derivative of Hermite polynomial of n'th degree
    if(m_dim == 0) {
        double sum = 0;
        for(int lambda=0; lambda<m_basisSize; lambda++) {
            sum += m_coefficients(lambda, n) * m_basis->evaluateSecondDerivative(x, lambda);
        }
        return sum;
    }
    else {
        return m_basis->evaluateSecondDerivative(x, n);
    }
}

double HermiteExpansion::basisElement(const int n, Eigen::VectorXd positions) {
    double prod = 1;
    for(int i=0; i<m_numberOfDimensions; i++) {
        m_dim = i;
        prod *= evaluate(positions(i), m_listOfStates(n, i));
    }
    return prod;
}

double HermiteExpansion::basisElementDer(const int n, const int i, Eigen::VectorXd positions) {
    // i is the dimension we are derivating with respect to
    double prod = evaluateDerivative(positions(i), m_listOfStates(n, i));
    for(int j=0; j<m_numberOfDimensions; j++) {
        if(i != j) {
            m_dim = j;
            prod *= evaluate(positions(j), m_listOfStates(n, j));
        }
    }
    return prod;
}

double HermiteExpansion::basisElementSecDer(const int n, const int i, Eigen::VectorXd positions) {
    // i is the dimension we are derivating with respect to
    double prod = evaluateSecondDerivative(positions(i), m_listOfStates(n, i));
    for(int j=0; j<m_numberOfDimensions; j++) {
        if(i != j) {
            m_dim = j;
            prod *= evaluate(positions(j), m_listOfStates(n, j));
        }
    }
    return prod;
}
