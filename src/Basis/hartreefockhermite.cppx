#include "hartreefockhermite.h"
#include "hermite.h"
#include "../system.h"
#include <iostream>
#include <fstream>

HartreeFockHermite::HartreeFockHermite(System *system, Basis *basis)  :
    Basis(system) {
    m_system                = system;
    m_numberOfParticles     = m_system->getNumberOfParticles();
    m_numberOfDimensions    = m_system->getNumberOfDimensions();
    m_omega                 = m_system->getFrequency();
    m_omegaSqrt             = sqrt(m_omega);
    m_path                  = m_system->getPath();
    m_basis                 = basis;
    readCoefficientFile();
    numberOfOrbitals();
    generateListOfStates(m_numberOfOrbitals);
}

std::ifstream::pos_type fileLength(std::string fileName)
{
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
    m_basisSize     = fileLength(fileName);
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

void HartreeFockHermite::numberOfOrbitals() {
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

void HartreeFockHermite::generateListOfStates(int orbitals) {
    m_basis->generateListOfStates(orbitals);
}

double HartreeFockHermite::basisElement(const int n, Eigen::VectorXd positions) {
    double sum = 0;
    for(int lambda=0; lambda<m_basisSize; lambda++) {
        sum += m_coefficients(n, lambda) * m_basis->basisElement(lambda, positions);
    }
    return sum;
}

double HartreeFockHermite::basisElementDer(const int n, const int i, Eigen::VectorXd positions) {
    // i is the dimension we are derivating with respect to
    double sum = 0;
    for(int lambda=0; lambda<m_basisSize; lambda++) {
        sum += m_coefficients(n, lambda) * m_basis->basisElementDer(lambda, i, positions);
    }
    return sum;
}

double HartreeFockHermite::basisElementSecDer(const int n, const int i, Eigen::VectorXd positions) {
    // i is the dimension we are derivating with respect to
    double sum = 0;
    for(int lambda=0; lambda<m_basisSize; lambda++) {
        sum += m_coefficients(n, lambda) * m_basis->basisElementSecDer(lambda, i, positions);
    }
    return sum;
}
