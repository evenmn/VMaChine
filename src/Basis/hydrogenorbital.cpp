#include "hydrogenorbital.h"
#include "../system.h"
#include <iostream>

HydrogenOrbital::HydrogenOrbital(System *system)  :
    Basis(system) {
    m_system                = system;
    m_Z                     = m_system->getAtomicNumber();
    m_numberOfParticles     = m_system->getNumberOfParticles();
    m_numberOfDimensions    = m_system->getNumberOfDimensions();
    assert(m_numberOfDimensions == 3);
    numberOfOrbitals();
}

double laguerre(const double x, const int n) {
    if(n == 0) {
        return 1;
    }
    else if(n == 1) {
        return 1 - x;
    }
    else {
        return ((2*n+1-x)*laguerre(x, n-1) - (n-1) * laguerre(x, n-2))/n;
    }
}

double associatedLaguerre(const double x, const unsigned int p, const unsigned int q) {
    if(q == 0) {
        return 1;
    }
    else if(q == 1) {
        return 1 + p - x;
    }
    else {
        return ((2*q+1+p-x)*associatedLaguerre(x, p, q-1) - (p+q-1) * associatedLaguerre(x, p, q-2))/q;
    }
}

unsigned int maxElectrons(const unsigned int i) {
    if(i==0) {
        return 2;
    }
    else {
        return 4 + maxElectrons(i-1);
    }
}

void HydrogenOrbital::numberOfOrbitals() {
    //Number of closed-shell orbitals
    int             orbital           = 0;
    unsigned int    i                 = 0;
    unsigned int    numberOfElectrons = 0;
    while(true) {
        for(unsigned int j=0; j<i; j++) {
            numberOfElectrons += maxElectrons(j);
            if(numberOfElectrons == m_numberOfParticles) {
                m_numberOfOrbitals = orbital+1;
                break;
            }
            else if(numberOfElectrons > m_numberOfParticles){
                std::cout << "This program supports closed-shells only. Please choose a P such that the orbital is full" << std::endl;
                exit(0);
            }
            orbital++;
        }
        i++;
    }
}

double HydrogenOrbital::evaluate(const double x, const unsigned int n) {
    //Hydrogen-like orbitals of a given n and l=0 (S-wave)
    double prefactor = (2*m_Z/n) * sqrt(2*m_Z/n) * sqrt(Basis::factorial(n-1)/(2 * n * Basis::factorial(n)));
    return prefactor * associatedLaguerre(2*m_Z*x/n, 1, n-1) * exp(-m_Z*x/n);
}

double HydrogenOrbital::evaluateDerivative(const double x, const unsigned int n) {
    //First derivative of Hydrogen-like orbitals of a given n and l=0 (S-wave)
    double prefactor = (2*m_Z/n) * sqrt(2*m_Z/n) * sqrt(Basis::factorial(n-1)/(2 * n * Basis::factorial(n)));
    return - prefactor * (associatedLaguerre(2*m_Z*x/n, 2, n-2)*exp(-m_Z*x/n) + associatedLaguerre(2*m_Z*x/n, 1, n-1)*exp(-m_Z*x/n) * m_Z/n);
}

double HydrogenOrbital::basisElement(const unsigned int n, const Eigen::VectorXd positions) {
    double r = positions.norm();
    double prefactor = (2*m_Z/n) * sqrt(2*m_Z/n) * sqrt(Basis::factorial(n-1)/(2 * n * Basis::factorial(n)));
    return prefactor * associatedLaguerre(2*m_Z*r/n, 1, n-1) * exp(-m_Z*r/n);
}

double HydrogenOrbital::basisElementDer(const unsigned int n, const unsigned int i, const Eigen::VectorXd positions) {
    // i is the dimension we are derivating with respect to
    double r = positions.norm();
    double prefactor = (2*m_Z/n) * sqrt(2*m_Z/n) * sqrt(Basis::factorial(n-1)/(2 * n * Basis::factorial(n)));
    return - prefactor * (associatedLaguerre(2*m_Z*r/n, 2, n-2)*exp(-m_Z*r/n) + associatedLaguerre(2*m_Z*r/n, 1, n-1)*exp(-m_Z*r/n) * m_Z/n);
}

double HydrogenOrbital::basisElementSecDer(const unsigned int n, const unsigned int i, const Eigen::VectorXd positions) {
    // i is the dimension we are derivating with respect to
    return 0;
}

void HydrogenOrbital::generateListOfStates(const int orbitals) {
    Eigen::MatrixXi listOfStates = Eigen::MatrixXi::Zero(m_numberOfParticles/2, m_numberOfDimensions);
    int counter = 0;
    // Three dimensions
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
