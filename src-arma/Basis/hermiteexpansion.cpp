#include "hermiteexpansion.h"
#include "../system.h"
#include "hermite.h"
#include <fstream>
#include <iostream>

HermiteExpansion::HermiteExpansion(System *system)
    : Basis(system)
{}

void HermiteExpansion::initialize()
{
    m_numberOfParticles = m_system->getNumberOfParticles();
    m_numberOfDimensions = m_system->getNumberOfDimensions();
    m_omega = m_system->getFrequency();
    m_path = m_system->getPath();
    m_basis = new Hermite(m_system);
    readCoefficientFile();
    Basis::numberOfOrbitals();
    Basis::generateListOfStates();

    arma::uword Ngrid = 1000;
    arma::vec x = arma::linspace(-10, 10, Ngrid);
    m_dim = 0;
    arma::vec y(Ngrid, arma::fill::zeros);
    for (arma::uword i = 0; i < Ngrid; i++) {
        y(i) = evaluate(x(i), 2) * exp(-0.5 * x(i) * x(i));
    }
    std::ofstream file;
    file.open(m_path + "test.dat");
    file << y << std::endl;
}

void HermiteExpansion::setParameters(arma::vec parameters) {}

std::string HermiteExpansion::generateFileName()
{
    std::string fileName = m_path;
    fileName += "arma::uword1/";
    fileName += "doublewell/";
    fileName += std::to_string(m_omega) + "w/";
    fileName += "coeffs.dat";
    return fileName;
}

void HermiteExpansion::readCoefficientFile()
{
    std::string fileName = generateFileName();
    m_basisSize = Basis::fileLength(fileName);
    m_coefficients.zeros(m_basisSize, m_basisSize);
    Basis::writeFileContent(fileName, m_coefficients);
}

double HermiteExpansion::evaluate(double x, arma::uword n)
{
    //Hermite polynomial of n'th degree
    if (m_dim == 0) {
        double sum = 0;
        for (arma::uword lambda = 0; lambda < m_basisSize; lambda++) {
            sum += m_coefficients(lambda, n) * m_basis->evaluate(x, lambda);
            //std::cout << m_coefficients(lambda, n) << std::endl;
        }
        return sum;
    } else {
        return m_basis->evaluate(x, n);
    }
}

double HermiteExpansion::evaluateDerivative(double x, arma::uword n)
{
    //First derivative of Hermite polynomial of n'th degree
    if (m_dim == 0) {
        double sum = 0;
        for (arma::uword lambda = 0; lambda < m_basisSize; lambda++) {
            sum += m_coefficients(lambda, n) * m_basis->evaluateDerivative(x, lambda);
        }
        return sum;
    } else {
        return m_basis->evaluateDerivative(x, n);
    }
}

double HermiteExpansion::evaluateSecondDerivative(const double x, const arma::uword n)
{
    //Second derivative of Hermite polynomial of n'th degree
    if (m_dim == 0) {
        double sum = 0;
        for (arma::uword lambda = 0; lambda < m_basisSize; lambda++) {
            sum += m_coefficients(lambda, n) * m_basis->evaluateSecondDerivative(x, lambda);
        }
        return sum;
    } else {
        return m_basis->evaluateSecondDerivative(x, n);
    }
}

double HermiteExpansion::basisElement(const arma::uword n, arma::vec positions)
{
    double prod = 1;
    for (m_dim = 0; m_dim < m_numberOfDimensions; m_dim++) {
        prod *= evaluate(positions(m_dim), m_listOfStates(n, m_dim));
    }
    return prod;
}

double HermiteExpansion::basisElementDer(const arma::uword n, const arma::uword i, arma::vec positions)
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

double HermiteExpansion::basisElementSecDer(const arma::uword n, const arma::uword i, arma::vec positions)
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

double HermiteExpansion::basisElementPar(const arma::uword n, arma::vec position)
{
    return 0;
}
