#include "hydrogenorbital.h"
#include "../system.h"

HydrogenOrbital::HydrogenOrbital(int Z)  :
    Basis() {
    m_Z = Z;
}

int fact(int n) {
    return (n == 1 || n == 0) ? 1 : fact(n - 1) * n;
}

double laguerre(double x, int n) {
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

double HydrogenOrbital::evaluate(double x, int n) {
    //Hydrogen-like orbitals of a given n and l=0 (S-wave)
    double prefactor = (2*m_Z/n) * sqrt(2*m_Z/n) * sqrt(fact(n-1)/(2 * n * fact(n)));
    return prefactor * associatedLaguerre(2*m_Z*x/n, 1, n-1) * exp(-m_Z*x/n);
}

double HydrogenOrbital::evaluateDerivative(double x, int n) {
    //First derivative of Hydrogen-like orbitals of a given n and l=0 (S-wave)
    double prefactor = (2*m_Z/n) * sqrt(2*m_Z/n) * sqrt(fact(n-1)/(2 * n * fact(n)));
    return - prefactor * (associatedLaguerre(2*m_Z*x/n, 2, n-2)*exp(-m_Z*x/n) + associatedLaguerre(2*m_Z*x/n, 1, n-1)*exp(-m_Z*x/n) * m_Z/n);
}
