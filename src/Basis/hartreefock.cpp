#include "hartreefock.h"
#include "../system.h"
#include "hermite.h"
#include "hermiteexpansion.h"
#include "hermitespin.h"
#include "hydrogenorbital.h"
#include "none.h"
#include <fstream>
#include <iostream>

HartreeFock::HartreeFock(System *system, Basis *basis)
    : Basis(system)
{
    m_system = system;
    m_numberOfParticles = m_system->getNumberOfParticles();
    m_numberOfDimensions = m_system->getNumberOfDimensions();
    m_omega = m_system->getFrequency();
    m_omegaSqrt = sqrt(m_omega);
    m_path = m_system->getPath();
    m_basis = basis;
    readCoefficientFile();
    Basis::numberOfOrbitals();
    Basis::generateListOfStates();
}

std::string HartreeFock::generateFileName()
{
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

void HartreeFock::readCoefficientFile()
{
    std::string fileName = generateFileName();
    std::ifstream inFile(generateFileName().c_str(), std::ios::in);
    m_basisSize = Basis::fileLength(fileName);
    m_coefficients = Eigen::MatrixXd::Zero(int(m_numberOfParticles / 2), m_basisSize);
    if (!inFile.is_open()) {
        std::cout << "file not found";
        MPI_Finalize();
        exit(0);
    } else {
        double value;
        int counter = 0;
        while (inFile >> value) {
            m_coefficients(int(counter / m_basisSize), counter % m_basisSize) = value;
            counter += 1;
        }
    }
}

void HartreeFock::setParameters(Eigen::VectorXd parameters) {}

double HartreeFock::evaluate(double x, int n)
{
    //Hermite polynomial of n'th degree
    double sum = 0;
    for (int lambda = 0; lambda < m_basisSize; lambda++) {
        sum += m_coefficients(n, lambda) * m_basis->evaluate(x, lambda);
    }
    return sum;
}

double HartreeFock::evaluateDerivative(double x, int n)
{
    //First derivative of Hermite polynomial of n'th degree
    double sum = 0;
    for (int lambda = 0; lambda < m_basisSize; lambda++) {
        sum += m_coefficients(n, lambda) * m_basis->evaluateDerivative(x, lambda);
    }
    return sum;
}

double HartreeFock::evaluateSecondDerivative(const double x, const int n)
{
    //Second derivative of Hermite polynomial of n'th degree
    double sum = 0;
    for (int lambda = 0; lambda < m_basisSize; lambda++) {
        sum += m_coefficients(n, lambda) * m_basis->evaluateSecondDerivative(x, lambda);
    }
    return sum;
}

double HartreeFock::basisElement(const int n, Eigen::VectorXd positions)
{
    double prod = 1;
    for (int i = 0; i < m_numberOfDimensions; i++) {
        prod *= evaluate(positions(i), int(m_listOfStates(n, i)));
    }
    return prod;
}

double HartreeFock::basisElementDer(const int n, const int i, Eigen::VectorXd positions)
{
    // i is the dimension we are derivating with respect to
    double prod = evaluateDerivative(positions(i), m_listOfStates(n, i));
    for (int j = 0; j < m_numberOfDimensions; j++) {
        if (i != j) {
            prod *= evaluate(positions(j), m_listOfStates(n, j));
        }
    }
    return prod;
}

double HartreeFock::basisElementSecDer(const int n, const int i, Eigen::VectorXd positions)
{
    // i is the dimension we are derivating with respect to
    double prod = evaluateSecondDerivative(positions(i), m_listOfStates(n, i));
    for (int j = 0; j < m_numberOfDimensions; j++) {
        if (i != j) {
            prod *= evaluate(positions(j), m_listOfStates(n, j));
        }
    }
    return prod;
}

double HartreeFock::basisElementPar(const int n, Eigen::VectorXd position)
{
    return 0;
}
