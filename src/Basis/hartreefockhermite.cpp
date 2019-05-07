#include "hartreefockhermite.h"
#include "hermite.h"
#include "../system.h"
#include <iostream>
#include <fstream>

HartreeFockHermite::HartreeFockHermite(System *system)  :
    Basis(system) {
    m_system                = system;
    m_numberOfParticles     = m_system->getNumberOfParticles();
    m_numberOfDimensions    = m_system->getNumberOfDimensions();
    m_omega                 = m_system->getFrequency();
    m_omegaSqrt             = sqrt(m_omega);
    m_path                  = m_system->getPath();
    m_hermite               = new Hermite(system);
    readCoefficientFile();
    numberOfOrbitals();
    generateListOfStates(m_numberOfOrbitals);
}

long fileLength(const std::string fileName) {
    std::ifstream inFile(fileName.c_str());
    return std::count(std::istreambuf_iterator<char>(inFile),
                      std::istreambuf_iterator<char>(), '\n');
}

std::string HartreeFockHermite::generateFileName() {
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

void HartreeFockHermite::readCoefficientFile() {
    std::string fileName = generateFileName();
    std::ifstream inFile(generateFileName().c_str(), std::ios::in);
    m_basisSize             = fileLength(fileName);
    m_coefficients          = Eigen::MatrixXd::Zero(int(m_numberOfParticles/2), m_basisSize);
    if (!inFile.is_open()) {
        std::cout << "file not found";
        MPI_Finalize();
        exit(0);
    }
    else {
        double value;
        unsigned int counter = 0;
        while (inFile >> value) {
            m_coefficients(int(counter/m_basisSize), counter % m_basisSize) = value;
            counter += 1;
        }
    }
}

void HartreeFockHermite::numberOfOrbitals() {
    int counter = 0;
    double orb = 0;
    while(orb <= m_basisSize * 2) {
        orb = 2 * Basis::binomial(counter, m_numberOfDimensions);
        m_numberOfOrbitals = counter+1;
        counter += 1;
    }
}

void HartreeFockHermite::generateListOfStates(const int orbitals) {
    m_hermite->generateListOfStates(orbitals);
}

double HartreeFockHermite::basisElement(const unsigned int n, const Eigen::VectorXd positions) {
    double sum = 0;
    for(unsigned int lambda=0; lambda<m_basisSize; lambda++) {
        sum += m_coefficients(n, lambda) * m_hermite->basisElement(lambda, positions);
    }
    return sum;
}

double HartreeFockHermite::basisElementDer(const unsigned int n, const unsigned int i, const Eigen::VectorXd positions) {
    // i is the dimension we are derivating with respect to
    double sum = 0;
    for(unsigned int lambda=0; lambda<m_basisSize; lambda++) {
        sum += m_coefficients(n, lambda) * m_hermite->basisElementDer(lambda, i, positions);
    }
    return sum;
}

double HartreeFockHermite::basisElementSecDer(const unsigned int n, const unsigned int i, const Eigen::VectorXd positions) {
    // i is the dimension we are derivating with respect to
    double sum = 0;
    for(unsigned int lambda=0; lambda<m_basisSize; lambda++) {
        sum += m_coefficients(n, lambda) * m_hermite->basisElementSecDer(lambda, i, positions);
    }
    return sum;
}
