#include "hermiteexpansion.h"
#include "../system.h"
#include "hermite.h"
#include <fstream>
#include <iostream>

HermiteExpansion::HermiteExpansion(System *system)
    : Basis(system)
{
    m_system = system;
    m_numberOfParticles = m_system->getNumberOfParticles();
    m_numberOfDimensions = m_system->getNumberOfDimensions();
    m_omega = m_system->getFrequency();
    m_path = m_system->getPath();
    m_basis = new Hermite(system);
    readCoefficientFile();
    Basis::numberOfOrbitals();
    Basis::generateListOfStates();

    int Ngrid = 1000;
    Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(Ngrid, -10, 10);
    m_dim = 0;
    Eigen::VectorXd y = Eigen::VectorXd::Zero(Ngrid);
    for (int i = 0; i < Ngrid; i++) {
        y(i) = evaluate(x(i), 2) * exp(-0.5 * x(i) * x(i));
    }
    std::ofstream file;
    file.open(m_path + "test.dat");
    file << y << std::endl;
}

void HermiteExpansion::setParameters(Eigen::VectorXd parameters) {}

std::string HermiteExpansion::generateFileName()
{
    std::string fileName = m_path;
    fileName += "int1/";
    fileName += "doublewell/";
    fileName += std::to_string(m_omega) + "w/";
    fileName += "coeffs.dat";
    return fileName;
}

void HermiteExpansion::readCoefficientFile()
{
    std::string fileName = generateFileName();
    m_basisSize = Basis::fileLength(fileName);
    m_coefficients = Eigen::MatrixXd::Zero(m_basisSize, m_basisSize);
    Basis::writeFileContentIntoEigenMatrix(fileName, m_coefficients);
}

double HermiteExpansion::evaluate(double x, int n)
{
    //Hermite polynomial of n'th degree
    if (m_dim == 0) {
        double sum = 0;
        for (int lambda = 0; lambda < m_basisSize; lambda++) {
            sum += m_coefficients(lambda, n) * m_basis->evaluate(x, lambda);
            //std::cout << m_coefficients(lambda, n) << std::endl;
        }
        return sum;
    } else {
        return m_basis->evaluate(x, n);
    }
}

double HermiteExpansion::evaluateDerivative(double x, int n)
{
    //First derivative of Hermite polynomial of n'th degree
    if (m_dim == 0) {
        double sum = 0;
        for (int lambda = 0; lambda < m_basisSize; lambda++) {
            sum += m_coefficients(lambda, n) * m_basis->evaluateDerivative(x, lambda);
        }
        return sum;
    } else {
        return m_basis->evaluateDerivative(x, n);
    }
}

double HermiteExpansion::evaluateSecondDerivative(const double x, const int n)
{
    //Second derivative of Hermite polynomial of n'th degree
    if (m_dim == 0) {
        double sum = 0;
        for (int lambda = 0; lambda < m_basisSize; lambda++) {
            sum += m_coefficients(lambda, n) * m_basis->evaluateSecondDerivative(x, lambda);
        }
        return sum;
    } else {
        return m_basis->evaluateSecondDerivative(x, n);
    }
}

double HermiteExpansion::basisElement(const int n, Eigen::VectorXd positions)
{
    double prod = 1;
    for (m_dim = 0; m_dim < m_numberOfDimensions; m_dim++) {
        prod *= evaluate(positions(m_dim), m_listOfStates(n, m_dim));
    }
    return prod;
}

double HermiteExpansion::basisElementDer(const int n, const int i, Eigen::VectorXd positions)
{
    // i is the dimension we are derivating with respect to
    double prod = evaluateDerivative(positions(i), m_listOfStates(n, i));
    for (m_dim = 0; m_dim < m_numberOfDimensions; m_dim++) {
        if (i != m_dim) {
            prod *= evaluate(positions(m_dim), m_listOfStates(n, m_dim));
        }
    }
    return prod;
}

double HermiteExpansion::basisElementSecDer(const int n, const int i, Eigen::VectorXd positions)
{
    // i is the dimension we are derivating with respect to
    double prod = evaluateSecondDerivative(positions(i), m_listOfStates(n, i));
    for (m_dim = 0; m_dim < m_numberOfDimensions; m_dim++) {
        if (i != m_dim) {
            prod *= evaluate(positions(m_dim), m_listOfStates(n, m_dim));
        }
    }
    return prod;
}

double HermiteExpansion::basisElementPar(const int n, Eigen::VectorXd position)
{
    return 0;
}
