//#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
#include "hermite.h"
#include "../system.h"
#include <iostream>

Hermite::Hermite(System *system)  :
    Basis(system) {
    m_system                = system;
    m_numberOfParticles     = m_system->getNumberOfParticles();
    m_numberOfDimensions    = m_system->getNumberOfDimensions();
    m_omega                 = m_system->getFrequency();
    m_omegaSqrt             = sqrt(m_omega);
    numberOfOrbitals();
}

int factorial(const int n) {
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

double binomial(const int n, const int p) {
    //Binomial coefficients, equal to magic numbers
    return factorial(n+p)/(factorial(n)*factorial(p));
}

void Hermite::numberOfOrbitals() {
    //Number of closed-shell orbitals
    int counter = 0;
    while(true) {
        double orb = 2*binomial(counter, m_numberOfDimensions);
        if(int(orb) == m_numberOfParticles) {
            m_numberOfOrbitals = counter+1;
            break;
        }
        else if(orb > m_numberOfParticles) {
            std::cout << "This program supports closed-shells only. Please choose a number of particles such that the orbital is full" << std::endl;
            exit(0);
        }
        counter += 1;
    }
}

double Hermite::evaluate(double x, int n) {
    //Hermite polynomial of n'th degree
    if(n == 0) {
        return 1;
    }
    else if(n == 1) {
        return 2 * m_omegaSqrt * x;
    }
    else {
        return 2 * (m_omegaSqrt * x * evaluate(x,n-1) - (n-1) * evaluate(x,n-2));
    }
    //return std::hermite(unsigned(n),x);
}

double Hermite::evaluateDerivative(double x, int n) {
    //First derivative of Hermite polynomial of n'th degree
    if(n == 0) {
        return 0;
    }
    else {
        return 2 * m_omegaSqrt * n * evaluate(x,n-1);
    }
}

double Hermite::evaluateSecondDerivative(double x, int n) {
    //Second derivative of Hermite polynomial of n'th degree
    if(n < 2) {
        return 0;
    }
    else {
        return 4 * m_omegaSqrt * n * (n-1) * evaluate(x,n-2);
    }
}

Eigen::MatrixXd Hermite::generateListOfStates() {
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
