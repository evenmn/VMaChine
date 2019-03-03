#include "hermite.h"
#include "../system.h"
#include <iostream>

Hermite::Hermite(System *system)  :
    Basis(system) {
    m_system                = system;
    m_numberOfParticles     = m_system->getNumberOfParticles();
    m_numberOfDimensions    = m_system->getNumberOfDimensions();
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
            std::cout << "This program supports closed-shells only. Please choose a P such that the orbital is full" << std::endl;
            exit(0);
        }
        counter += 1;
    }
}

double Hermite::evaluate(double x, int n) {
    //Hermite polynomial of n'th degree
    //std::cout << n << std::endl;
    if(n == 0) {
        return 1;
    }
    else if(n == 1) {
        return 2*x;
    }
    else {
        return 2*x*evaluate(x,n-1)-2*(n-1)*evaluate(x,n-2);
    }
}

double Hermite::evaluateDerivative(double x, int n) {
    //First derivative of Hermite polynomial of n'th degree
    if(n == 0) {
        return 0;
    }
    else {
        return 2*n*evaluate(x,n-1);
    }
}
