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
    Basis::generateListOfStates();
}

double associatedLaguerre(double x, int p, int q) {
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

double associatedLaguerreDer(double x, int p, int q, int k) {
    if(k <= q) {
        return pow(-1,k)*associatedLaguerre(x, p+k, q-k);
    }
    else {
        return 0;
    }
}

double HydrogenOrbital::radial(double r, int n) {
    double prefactor = 1;// (2*m_Z/n) * sqrt(2*m_Z/n) * sqrt(Basis::factorial(n-1)/(2 * n * Basis::factorial(n)));
    return prefactor * associatedLaguerre(2*m_Z*r/n, 1, n-1) * exp(-m_Z*r/n);
}

double HydrogenOrbital::angular(double theta, double phi, int l, int m) {
    return theta + phi + l + m;
}

int maxElectrons(int i) {
    if(i==0) {
        return 2;
    }
    else {
        return 4 + maxElectrons(i-1);
    }
}

double HydrogenOrbital::evaluate(double r, int n) {
    //Hydrogen-like orbitals of a given n and l=0 (S-wave)
    return radial(r, n);
}

double HydrogenOrbital::evaluateDerivative(double r, int n) {
    //First derivative of Hydrogen-like orbitals of a given n and l=0 (the S-wave)
    double prefactor = -m_Z/n;
    return prefactor * radial(r, n) + exp(-m_Z*r/n)*associatedLaguerreDer(2*m_Z*r/n, 1, n-1, 1);
}

double HydrogenOrbital::evaluateSecondDerivative(double r, int n) {
    //Second derivative of Hydrogen-like orbitals of a given n and l=0 (the S-wave)
    double prefactor = m_Z * m_Z/(n * n);
    return prefactor * radial(r, n) + exp(-m_Z*r/n)*associatedLaguerreDer(2*m_Z*r/n, 1, n-1, 2);
}

double HydrogenOrbital::basisElement(const int n, Eigen::VectorXd positions) {
    return evaluate(positions.norm(), n+1);
}

double HydrogenOrbital::basisElementDer(const int n, const int i, Eigen::VectorXd positions) {
    double r = positions.norm();
    return (positions(i)/r)*evaluateDerivative(r, n+1);
}

double HydrogenOrbital::basisElementSecDer(const int n, const int i, Eigen::VectorXd positions) {
    double r = positions.norm();
    return (1/r-(positions(i)*positions(i))/(r*r*r))*evaluateSecondDerivative(r,n+1);
}
