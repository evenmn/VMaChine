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
    numberOfOrbitalss();
    generateLOS();
}

void HydrogenOrbital::setParameters(Eigen::VectorXd parameters) {
    m_alpha = parameters(0);
}

void HydrogenOrbital::numberOfOrbitalss() {
    int number = 1;
    int number2 = 0;
    int maxNumber = 0;
    int counter = 1;
    m_numberOfShells = 1;
    while(true) {
        number2 += 2*number;
        if(m_numberOfParticles==number2) {
            m_numberOfOrbitalss = counter;
            break;
        }
        else if(m_numberOfParticles<number2) {
            std::cout << "Only closed shell atoms are accepted: He, Be, Ne, Mg, Ar, Zn, Kr" << std::endl;
            MPI_Finalize();
            exit(0);
        }
        else {
            if(number > maxNumber) {
                m_numberOfShells += 1;
                maxNumber = number;
                number = 1;
            }
            else {
                number += 2;
            }
        }
        counter++;
    }
}

void HydrogenOrbital::generateLOS() {
    int sum = 0;
    for(int i=1; i<m_numberOfShells+1; i++) {
        sum += i*i;
    }
    m_LOS = Eigen::MatrixX3i::Zero(sum,3);
    int i=0;
    for(int n=1; n<m_numberOfShells+1; n++) {
        for(int l=0; l<n; l++) {
            for(int m=-l; m<l+1; m++) {
                m_LOS.row(i) << n, l, m;
                i++;
            }
        }
    }
}

double HydrogenOrbital::basisElement(const int n, Eigen::VectorXd position) {
    return evaluateCart(position, m_LOS(n,0), m_LOS(n,1), m_LOS(n,2));
}

double HydrogenOrbital::basisElementDer(const int n, const int i, Eigen::VectorXd position) {
    return evaluateCartDerivative(position, i, m_LOS(n,0), m_LOS(n,1), m_LOS(n,2));
}

double HydrogenOrbital::basisElementSecDer(const int n, const int i, Eigen::VectorXd position) {
    return evaluateCartSecondDerivative(position, i, m_LOS(n,0), m_LOS(n,1), m_LOS(n,2));
}

double HydrogenOrbital::basisElementPar(const int n, Eigen::VectorXd position) {
    double r = position.norm();
    if(n==0) {
        return -m_Z*r;
    }
    else if(n==1) {
        return m_Z*r*(2-m_alpha*m_Z*r/2);
    }
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
            result = 2-k*r;
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
    return result * exp(-k*r/n);
}

double HydrogenOrbital::evaluateCartDerivative(Eigen::VectorXd position, int i, int n, int l, int m) {
    //Hard coded derivatives of Hydrogen orbitals taken from Jorgen Hogberget
    double r = position.norm();
    double x = position(0);
    double y = position(1);
    double z = position(2);
    double t = position(i);
    double k = m_alpha*m_Z;
    double result = 0;
    if(n==1) {
        result = -k*t;
    }
    else if(n==2) {
        if(l==0) {
            result = -k*t*(4-k*r);
        }
        else if(l==1) {
            if(m==0) {
                result = -k*t*z;
                if(i==2) {
                    result += 2*r;
                }
            }
            else if(m==1) {
                result = -k*t*x;
                if(i==0) {
                    result += 2*r;
                }
            }
            else if(m==-1) {
                result = -k*t*y;
                if(i==1) {
                    result += 2*r;
                }
            }
        }
    }
    else if(n==3) {
        if(l==0) {
            result = -k*t*(2*k*k*r*r - 30*k*r + 81);
        }
        else if(l==1) {
            if(m==0) {
                if(i==0) {
                    result = -k*t*z*(k*r-9);
                }
                else if(i==1) {
                    result = -k*t*z*(k*r-9);
                }
                else if(i==2) {
                    result = -k*t*z*(k*r-9) - 18*r + 3*k*r*r;
                }
            }
            else if(m==1) {
                if(i==0) {
                    result = -k*t*x*(k*r-9) - 18*r + 3*k*r*r;
                }
                else if(i==1) {
                    result = -k*t*x*(k*r-9);
                }
                else if(i==2) {
                    result = -k*t*x*(k*r-9);
                }
            }
            else if(m==-1) {
                if(i==0) {
                    result = -k*t*y*(k*r-9);
                }
                else if(i==1) {
                    result = -k*t*y*(k*r-9) - 18*r + 3*k*r*r;
                }
                else if(i==2) {
                    result = -k*t*y*(k*r-9);
                }
            }
        }
        else if(l==2) {
            if(m==0) {
                if(i==0) {
                    result = -t*(k*(-r*r+3*z*z)+6*r);
                }
                else if(i==1) {
                    result = -t*(k*(-r*r+3*z*z)+6*r);
                }
                else if(i==2) {
                    result = -t*(k*(-r*r+3*z*z)-12*r);
                }
            }
            else if(m==1) {
                if(i==0) {
                    result = -z*(k*x*x-3*r);
                }
                else if(i==1) {
                    result = -y*(k*x*y);
                }
                else if(i==2) {
                    result = -x*(k*z*z-3*r);
                }
            }
            else if(m==-1) {
                if(i==0) {
                    result = -k*x*y*z;
                }
                else if(i==1) {
                    result = -z*(k*y*y-3*r);
                }
                else if(i==2) {
                    result = -y*(k*z*z-3*r);
                }
            }
            else if(m==2) {
                if(i==0) {
                    result = -x*(k*(x*x-y*y)-6*r);
                }
                else if(i==1) {
                    result = -y*(k*(x*x-y*y)+6*r);
                }
                else if(i==2) {
                    result = -k*z*(x*x-y*y);
                }
            }
            else if(m==-2) {
                if(i==0) {
                    result = -y*(k*x*x-3*r);
                }
                else if(i==1) {
                    result = -x*(k*y*y-3*r);
                }
                else if(i==2) {
                    result = -k*x*y*z;
                }
            }
        }
    }
    else if(n==4) {
        if(l==0) {
            result = -k*t*(k*k*k*r*r*r-36*k*k*r*r+336*k*r-768);
        }
        else if(l==1) {
            if(m==0) {
                if(i==0) {
                    result = -k*t*z*(k*r-20)*(k*r-8);
                }
                else if(i==1) {
                    result = -k*t*z*(k*r-20)*(k*r-8);
                }
                else if(i==2) {
                    result = -k*t*z*(k*r-20)*(k*r-8)+320*r-80*k*r*r+4*k*k*r*r*r;
                }
            }
            else if(m==1) {
                if(i==0) {
                    result = -k*t*x*(k*r-20)*(k*r-8)+320*r-80*k*r*r+4*k*k*r*r*r;
                }
                else if(i==1) {
                    result = -k*t*x*(k*r-20)*(k*r-8);
                }
                else if(i==2) {
                    result = -k*t*x*(k*r-20)*(k*r-8);
                }
            }
            else if(m==-1) {
                if(i==0) {
                    result = -k*t*y*(k*r-20)*(k*r-8);
                }
                else if(i==1) {
                    result = -k*t*y*(k*r-20)*(k*r-8)+320*r-80*k*r*r+4*k*k*r*r*r;
                }
                else if(i==2) {
                    result = -k*t*y*(k*r-20)*(k*r-8);
                }
            }
        }
    }
    return result * exp(-k*r/n)/(n*r);
}

double HydrogenOrbital::evaluateCartSecondDerivative(Eigen::VectorXd position, int i, int n, int l, int m) {
    //Hard coded Laplacian of Hydrogen orbitals
    double r = position.norm();
    double x = position(0);
    double y = position(1);
    double z = position(2);
    double t = position(i);
    double k = m_alpha*m_Z;

    double result = 0;
    if(n==1) {
        result = k*(t*t*(k*r+1)-r*r);
    }
    else if(n==2) {
        if(l==0) {
            result = -k*(t*t*(k*r*(k*r-4)-8)-2*r*r*(k*r-4));
        }
        else if(l==1) {
            if(m==-1) {
                result = -k*y*(-t*t*(k*r+2) + 2*r*r);
                if(i==1) {
                    result -= 4*k*y*r*r;
                }
            }
            else if(m==0) {
                result = -k*z*(-t*t*(k*r+2) + 2*r*r);
                if(i==2) {
                    result -= 4*k*z*r*r;
                }
            }
            else if(m==1) {
                result = -k*x*(-t*t*(k*r+2) + 2*r*r);
                if(i==0) {
                    result -= 4*k*x*r*r;
                }
            }
        }
    }
    else if(n==3) {
        if(l==0) {
            result = k*(r*r*(2*k*r*(5-3*k*r)-27)+t*t*(k*r*(2*k*r*(k*r-18)+9)+27))/9;
        }
        else if(l==1) {
            if(m==0) {
                result = k*z*(k*r-18)*(k*r-6)/(9*r);
            }
            else if(m==1) {
                result = k*x*(k*r-18)*(k*r-6)/(9*r);
            }
            else if(m==-1) {
                result = k*y*(k*r-18)*(k*r-6)/(9*r);
            }
        }
        else if(l==2) {
            if(m==0) {
                result = k*(-r*r+3*z*z)*(k*r-18)/(9*r);
            }
            else if(m==1) {
                result = k*x*z*(k*r-18)/(9*r);
            }
            else if(m==-1) {
                result = k*y*z*(k*r-18)/(9*r);
            }
            else if(m==2) {
                result = k*(x*x-y*y)*(k*r-18)/(9*r);
            }
            else if(m==-2) {
                result = k*x*y*(k*r-18)/(9*r);
            }
        }
    }
    else if(n==4) {
        if(l==0) {
            result = k*(k*r-32)*(k*k*k*r*r*r-24*k*k*r*r+144*k*r-192)/(16*r);
        }
        else if(l==1) {
            if(m==0) {
                result = k*z*(k*r-32)*(k*k*r*r-20*k*r+80)/(16*r);
            }
            else if(m==1) {
                result = k*x*(k*r-32)*(k*k*r*r-20*k*r+80)/(16*r);
            }
            else if(m==-1) {
                result = k*y*(k*r-32)*(k*k*r*r-20*k*r+80)/(16*r);
            }
        }
    }
    return result * exp(-k*r/n) / (n*n*r*r*r);
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
