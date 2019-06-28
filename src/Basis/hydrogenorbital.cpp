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

void HydrogenOrbital::numberOfOrbitalss() {
    int number = 1;
    int number2 = 2;
    int maxNumber = 0;
    int counter = 1;
    m_numberOfShells = 1;
    while(true) {
        std::cout << number2 << std::endl;
        if(m_numberOfParticles==number2) {
            m_numberOfOrbitalss = counter;
            break;
        }
        else if(m_numberOfParticles<2*number2) {
            std::cout << "What u tryna do?" << std::endl;
            MPI_Finalize();
            exit(0);
        }
        else {
            number2 += 2*number;
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
    m_LOS = Eigen::MatrixX3i::Zero(m_numberOfParticles/2,3);
    int i=0;
    for(int n=1; n<m_numberOfShells+1; n++) {
        for(int l=0; i<n; l++) {
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
    double k = m_alpha*m_Z;
    double result = 0;
    if(n==1) {
        result = -k*position(i)/r;
    }
    else if(n==2) {
        if(l==0) {
            result = -k*position(i)*(4-k*r)/(2*r);
        }
        else if(l==1) {
            if(m==0) {
                if(i==0) {
                    result = -k*x*z/(2*r);
                }
                else if(i==1) {
                    result = -k*y*z/(2*r);
                }
                else if(i==2) {
                    result = -k*z*z/(2*r) + 1;
                }
            }
            else if(m==1) {
                if(i==0) {
                    result = -k*x*x/(2*r) + 1;
                }
                else if(i==1) {
                    result = -k*y*x/(2*r);
                }
                else if(i==2) {
                    result = -k*z*x/(2*r);
                }
            }
            else if(m==-1) {
                if(i==0) {
                    result = -k*x*y/(2*r);
                }
                else if(i==1) {
                    result = -k*y*y/(2*r) + 1;
                }
                else if(i==2) {
                    result = -k*z*y/(2*r);
                }
            }
        }
    }
    else if(n==3) {
        if(l==0) {
            result = -k*position(i)*(2*k*k*r*r - 30*k*r + 81)/(3*r);
        }
        else if(l==1) {
            if(m==0) {
                if(i==0) {
                    result = -k*x*z*(k*r-9)/(3*r);
                }
                else if(i==1) {
                    result = -k*y*z*(k*r-9)/(3*r);
                }
                else if(i==2) {
                    result = -k*z*z*(k*r-9)/(3*r) - 6 + k*r;
                }
            }
            else if(m==1) {
                if(i==0) {
                    result = -k*x*z*(k*r-9)/(3*r) - 6 + k*r;
                }
                else if(i==1) {
                    result = -k*y*z*(k*r-9)/(3*r);
                }
                else if(i==2) {
                    result = -k*z*z*(k*r-9)/(3*r);
                }
            }
            else if(m==-1) {
                if(i==0) {
                    result = -k*x*z*(k*r-9)/(3*r);
                }
                else if(i==1) {
                    result = -k*y*z*(k*r-9)/(3*r) - 6 + k*r;
                }
                else if(i==2) {
                    result = -k*z*z*(k*r-9)/(3*r);
                }
            }
        }
        else if(l==2) {
            if(m==0) {
                if(i==0) {
                    result = -x*(k*(-r*r+3*z*z)+6*r)/(3*r);
                }
                else if(i==1) {
                    result = -y*(k*(-r*r+3*z*z)+6*r)/(3*r);
                }
                else if(i==2) {
                    result = -z*(k*(-r*r+3*z*z)-12*r)/(3*r);
                }
            }
            else if(m==1) {
                if(i==0) {
                    result = -z*(k*x*x-3*r)/(3*r);
                }
                else if(i==1) {
                    result = -k*x*y*z/(3*r);
                }
                else if(i==2) {
                    result = -x*(k*z*z-3*r)/(3*r);
                }
            }
            else if(m==-1) {
                if(i==0) {
                    result = -k*x*y*z/(3*r);
                }
                else if(i==1) {
                    result = -z*(k*y*y-3*r)/(3*r);
                }
                else if(i==2) {
                    result = -y*(k*z*z-3*r)/(3*r);
                }
            }
            else if(m==2) {
                if(i==0) {
                    result = -x*(k*(x*x-y*y)-6*r)/(3*r);
                }
                else if(i==1) {
                    result = -y*(k*(x*x-y*y)+6*r)/(3*r);
                }
                else if(i==2) {
                    result = -k*z*(x*x-y*y)/(3*r);
                }
            }
            else if(m==-2) {
                if(i==0) {
                    result = -y*(k*x*x-3*r)/(3*r);
                }
                else if(i==1) {
                    result = -x*(k*y*y-3*r)/(3*r);
                }
                else if(i==2) {
                    result = -k*x*y*z/(3*r);
                }
            }
        }
    }
    else if(n==4) {
        if(l==0) {
            result = -k*position(i)*(k*k*k*r*r*r-36*k*k*r*r+336*k*r-768)/(4*r);
        }
        else if(l==1) {
            if(m==0) {
                if(i==0) {
                    result = -k*x*z*(k*r-20)*(k*r-8)/(4*r);
                }
                else if(i==1) {
                    result = -k*y*z*(k*r-20)*(k*r-8)/(4*r);
                }
                else if(i==2) {
                    result = -k*z*z*(k*r-20)*(k*r-8)/(4*r)+320-80*k*r+4*k*k*r*r;
                }
            }
            else if(m==1) {
                if(i==0) {
                    result = -k*x*x*(k*r-20)*(k*r-8)/(4*r)+320-80*k*r+4*k*k*r*r;
                }
                else if(i==1) {
                    result = -k*x*y*(k*r-20)*(k*r-8)/(4*r);
                }
                else if(i==2) {
                    result = -k*x*z*(k*r-20)*(k*r-8)/(4*r);
                }
            }
            else if(m==-1) {
                if(i==0) {
                    result = -k*x*y*(k*r-20)*(k*r-8)/(4*r);
                }
                else if(i==1) {
                    result = -k*y*y*(k*r-20)*(k*r-8)/(4*r)+320-80*k*r+4*k*k*r*r;
                }
                else if(i==2) {
                    result = -k*y*z*(k*r-20)*(k*r-8)/(4*r);
                }
            }
        }
    }
    return result * exp(-k*r/n);
}

double HydrogenOrbital::evaluateCartSecondDerivative(Eigen::VectorXd position, int i, int n, int l, int m) {
    //Hard coded Laplacian of Hydrogen orbitals taken from Jorgen Hogberget
    double r = position.norm();
    double x = position(0);
    double y = position(1);
    double z = position(2);
    double k = m_alpha*m_Z;
    double result = 0;
    if(n==1) {
        if(i==0){
            result = k*(x*x*(k*r+1)-r*r)/(r*r*r);
        }
        else if(i==1){
            result = k*(y*y*(k*r+1)-r*r)/(r*r*r);
        }
        else if(i==2){
            result = k*(z*z*(k*r+1)-r*r)/(r*r*r);
        }
    }
    else if(n==2) {
        if(l==0) {
            if(i==0) {
                result = -k*(k*r*(x*x*(k*r-4)-2*r*r)+8*(y*y+z*z))/(4*r*r*r);
            }
            else if(i==1) {
                result = -k*(k*r*(y*y*(k*r-4)-2*r*r)+8*(x*x+z*z))/(4*r*r*r);
            }
            else if(i==2) {
                result = -k*(k*r*(z*z*(k*r-4)-2*r*r)+8*(x*x+y*y))/(4*r*r*r);
            }
        }
        else if(l==1) {
            if(m==-1) {
                if(i==0) {
                    result = k*y*(-k*r*x*x+2*y*y+2*z*z)/(4*r*r*r);
                }
                else if(i==1) {
                    result = k*y*(-k*r*y*y+4*y*y+6*x*x+6*z*z)/(4*r*r*r);
                }
                else if(i==2) {
                    result = k*y*(-k*r*z*z+2*x*x+2*y*y)/(4*r*r*r);
                }
            }
            else if(m==0) {
                if(i==0) {
                    result = k*z*(-k*r*x*x+2*y*y+2*z*z)/(4*r*r*r);
                }
                else if(i==1) {
                    result = k*z*(-k*r*y*y+2*x*x+2*z*z)/(4*r*r*r);
                }
                else if(i==2) {
                    result = k*z*(-k*r*z*z+4*z*z+6*x*x+6*y*y)/(4*r*r*r);
                }
            }
            else if(m==1) {
                if(i==0) {
                    result = k*x*(-k*r*x*x+4*x*x+6*y*y+6*z*z)/(4*r*r*r);
                }
                else if(i==1) {
                    result = k*x*(-k*r*y*y+2*x*x+2*z*z)/(4*r*r*r);
                }
                else if(i==2) {
                    result = k*x*(-k*r*z*z+2*x*x+2*y*y)/(4*r*r*r);
                }
            }
        }
    }
    else if(n==3) {
        if(l==0) {
            result = k*(k*r-18)*(2*k*k*r*r-18*k*r+27)/(9*r);
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
    return result * exp(-k*r/n);
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
