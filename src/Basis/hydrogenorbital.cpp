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

int maxElectrons(int i) {
    if(i==0) {
        return 2;
    }
    else {
        return 4 + maxElectrons(i-1);
    }
}

double HydrogenOrbital::basisElement(const int n, Eigen::VectorXd position) {
    double p = m_Z*position.norm()/(n+1);
    return exp(-p) * associatedLaguerre(2*p, 1, n);
    //return evaluateCart(position, n+1, 0, 0);
}

double HydrogenOrbital::basisElementDer(const int n, const int i, Eigen::VectorXd position) {
    double r = position.norm();
    double p = m_Z*r/(n+1);
    double q = m_Z * position(i) / ((n+1) * r);
    //return exp(-p) * (associatedLaguerreDer(2*p, 1, n, 1) - 2*p * (position(i)/(r*r)) * associatedLaguerre(2*p, 1, n));
    return q * exp(-p) * (2 * associatedLaguerreDer(2*p, 1, n, 1) - associatedLaguerre(2*p, 1, n));
}

double HydrogenOrbital::basisElementSecDer(const int n, const int i, Eigen::VectorXd position) {
    double r = position.norm();
    double p = m_Z*r/(n+1);
    double q = m_Z / ((n+1) * r);
    return q * exp(-p) * (2 * associatedLaguerreDer(2*p, 1, n, 1) * (1 - 2*q*position(i)*position(i) - position(i)*position(i)/(r*r)) + associatedLaguerre(2*p, 1, n)*(position(i)*position(i)/(r*r)+m_Z*position(i)*position(i)/((n+1)*r)-1) + 4 * q * position(i) * position(i) * associatedLaguerreDer(2*p, 1, n, 2));
    //double q = position(i)*position(i)/(r*r);
    //double s = 1 - q*(1+p);
    //return exp(-p) * (s*(2*associatedLaguerreDer(2*p,1,n,1) - associatedLaguerre(2*p,1,n)) + 4*p*q*associatedLaguerreDer(2*p,1,n,2))*m_Z/((n+1)*r);
}

double HydrogenOrbital::evaluateCart(Eigen::VectorXd position, int n, int l, int m) {
    //Hard coded Hydrogen orbitals taken from Jorgen Hogberget
    double r = position.norm();
    double x = position(0);
    double y = position(1);
    double z = position(2);
    double k = m_alpha*m_Z;
    double result = 0;
    if(n==1) {
        result = 1;
    }
    else if(n==2) {
        if(l==0) {
            result = k*r - 2;
        }
        else if(l==1) {
            if(m==0) {
                result = z;
            }
            else if(m==1) {
                result = x;
            }
            else if(m==-1) {
                result = y;
            }
        }
    }
    else if(n==3) {
        if(l==0) {
            result = 2*k*k*r*r-18*k*r+27;
        }
        else if(l==1) {
            if(m==0) {
                result = z*(k*r-6);
            }
            else if(m==1) {
                result = x*(k*r-6);
            }
            else if(m==-1) {
                result = y*(k*r-6);
            }
        }
        else if(l==2) {
            if(m==0) {
                result = -r*r+3*z*z;
            }
            else if(m==1) {
                result = x*z;
            }
            else if(m==-1) {
                result = y*z;
            }
            else if(m==2) {
                result = x*x - y*y;
            }
            else if(m==-2) {
                result = x*y;
            }
        }
    }
    else if(n==4) {
        if(l==0) {
            result = k*k*k*r*r*r-24*k*k*r*r+144*k*r-192;
        }
        else if(l==1) {
            if(m==0) {
                result = z*(k*k*r*r-20*k*r+80);
            }
            else if(m==1) {
                result = x*(k*k*r*r-20*k*r+80);
            }
            else if(m==-1) {
                result = y*(k*k*r*r-20*k*r+80);
            }
        }
    }
    return result * exp(-k*r);
}

double HydrogenOrbital::evaluateCartDerivative(Eigen::VectorXd position, int i, int n, int l, int m) {
    //First derivative of Hydrogen-like orbitals of a given n and l=0 (the S-wave)
    double r = position.norm();
    return r+n; //exp(-m_Z*r/n) * (associatedLaguerreDer(2*m_Z*r/n, 1, n-1, 1) - (m_Z/n) * associatedLaguerre(2*m_Z*r/n, 1, n-1));
}

double HydrogenOrbital::evaluateCartSecondDerivative(Eigen::VectorXd position, int i, int n, int l, int m) {
    //Second derivative of Hydrogen-like orbitals of a given n and l=0 (the S-wave)
    double r = position.norm();
    //double prefactor = m_Z/n;
    return position(0)+n; //prefactor * radial(r, n) + exp(-m_Z*r/n)*associatedLaguerreDer(2*m_Z*r/n, 1, n-1, 2);
}

double HydrogenOrbital::evaluate(double x, int n) {
    return x+n;
}

double HydrogenOrbital::evaluateDerivative(double x, int n) {
    return x+n;
}

double HydrogenOrbital::evaluateSecondDerivative(double x, int n) {
    return x+n;
}
