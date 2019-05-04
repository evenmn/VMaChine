#include "hartreefockhermite.h"
#include "hermite.h"
#include "../system.h"
#include <iostream>

HartreeFockHermite::HartreeFockHermite(System *system)  :
    Basis(system) {
    m_system                = system;
    m_numberOfParticles     = m_system->getNumberOfParticles();
    m_numberOfDimensions    = m_system->getNumberOfDimensions();
    m_omega                 = m_system->getFrequency();
    m_omegaSqrt             = sqrt(m_omega);
    m_basisSize             = 120;
    m_coefficients          = Eigen::MatrixXd::Zero(int(m_numberOfParticles/2), m_basisSize);
    m_hermite               = new Hermite(system);
    numberOfOrbitals();

}

void HartreeFockHermite::numberOfOrbitals() {
    //Number of closed-shell orbitals
    int counter = 0;
    while(true) {
        double orb = 2 * Basis::binomial(counter, m_numberOfDimensions);
        if(int(orb) == m_numberOfParticles) {
            m_numberOfOrbitals = counter+1;
            break;
        }
        else if(orb > m_numberOfParticles) {
            std::cout << "This program supports closed-shells only. Please choose a number of particles such that the orbital is full" << std::endl;
            MPI_Finalize();
            exit(0);
        }
        counter += 1;
    }
}

double HartreeFockHermite::evaluate(double x, int n) {
    //Hermite polynomial of n'th degree
    double sum = 0;
    for(unsigned int lambda=0; lambda<m_basisSize; lambda++) {
        sum += m_coefficients(n, lambda) * m_hermite->evaluate(x, lambda);
    }
    return sum;
}

double HartreeFockHermite::evaluateDerivative(double x, int n) {
    //First derivative of Hermite polynomial of n'th degree
    double sum = 0;
    for(unsigned int lambda=0; lambda<m_basisSize; lambda++) {
        sum += m_coefficients(n, lambda) * m_hermite->evaluateDerivative(x, lambda);
    }
    return sum;
}

double HartreeFockHermite::evaluateSecondDerivative(double x, int n) {
    //Second derivative of Hermite polynomial of n'th degree
    double sum = 0;
    for(unsigned int lambda=0; lambda<m_basisSize; lambda++) {
        sum += m_coefficients(n, lambda) * m_hermite->evaluateSecondDerivative(x, lambda);
    }
    return sum;
}

Eigen::MatrixXd HartreeFockHermite::generateListOfStates() {
    // Returns the index list used in Slater
    // For instance (0,0), (1,0), (0,1) for 6P in 2D
    //              (0,0,0), (1,0,0), (0,1,0), (0,0,1) for 8P in 3D etc..
    numberOfOrbitals();
    Eigen::MatrixXd listOfStates = Eigen::MatrixXd::Zero(m_numberOfParticles/2, m_numberOfDimensions);
    int counter = 0;
    // One dimension
    if (m_numberOfDimensions == 1) {
        for(int i=0; i<m_numberOfOrbitals; i++) {
            listOfStates(i) = i;
        }
    }
    // Two dimensions
    if (m_numberOfDimensions == 2) {
        for(int i=0; i<m_numberOfOrbitals; i++) {
            for(int s=i; s<m_numberOfOrbitals; s++) {
                int j = s - i;
                listOfStates(counter,1) = i;
                listOfStates(counter,0) = j;
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
                    listOfStates(counter,0) = i;
                    listOfStates(counter,1) = j;
                    listOfStates(counter,2) = k;
                    counter += 1;
                }
            }
        }
    }
    // Four dimensions
    else if (m_numberOfDimensions == 4) {
        for(int i=0; i<m_numberOfOrbitals; i++) {
            for(int j=0; j<m_numberOfOrbitals; j++) {
                for(int k=0; k<m_numberOfOrbitals; k++) {
                    for(int s=i+j; s<m_numberOfOrbitals; s++) {
                        int l = s - i - j - k;
                        listOfStates(counter,0) = i;
                        listOfStates(counter,1) = j;
                        listOfStates(counter,2) = k;
                        listOfStates(counter,3) = l;
                        counter += 1;
                    }
                }
            }
        }
    }

    else {
        std::cout << "Number of dimensions should be in the range [1, 4]" << std::endl;
        exit(0);
    }
    return listOfStates;
}
