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
    generateListOfStates(m_numberOfOrbitals);
}

void Hermite::numberOfOrbitals() {
    //Number of closed-shell orbitals
    int counter = 0;
    while(true) {
        unsigned int orb = 2 * Basis::binomial(unsigned(counter), m_numberOfDimensions);
        if(orb == m_numberOfParticles) {
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

void Hermite::generateListOfStates(const int orbitals) {
    // Returns the index list used in Slater
    // For instance (0,0), (1,0), (0,1) for 6P in 2D
    //              (0,0,0), (1,0,0), (0,1,0), (0,0,1) for 8P in 3D etc..
    unsigned int numberOfStates = Basis::binomial(unsigned(orbitals-1), unsigned(m_numberOfDimensions));
    m_listOfStates = Eigen::MatrixXi::Zero(numberOfStates, m_numberOfDimensions);
    unsigned int counter = 0;
    // One dimension
    if (m_numberOfDimensions == 1) {
        for(int i=0; i<orbitals; i++) {
            m_listOfStates(i) = i;
        }
    }
    // Two dimensions
    if (m_numberOfDimensions == 2) {
        for(int i=0; i<orbitals; i++) {
            for(int s=i; s<orbitals; s++) {
                int j = s - i;
                m_listOfStates(counter,1) = i;
                m_listOfStates(counter,0) = j;
                counter += 1;
            }
        }
    }
    // Three dimensions
    else if (m_numberOfDimensions == 3) {
        for(int i=0; i<orbitals; i++) {
            for(int j=0; j<orbitals; j++) {
                for(int s=i+j; s<orbitals; s++) {
                    int k = s - i - j;
                    m_listOfStates(counter,0) = i;
                    m_listOfStates(counter,1) = j;
                    m_listOfStates(counter,2) = k;
                    counter += 1;
                }
            }
        }
    }
    // Four dimensions
    else if (m_numberOfDimensions == 4) {
        for(int i=0; i<orbitals; i++) {
            for(int j=0; j<orbitals; j++) {
                for(int k=0; k<orbitals; k++) {
                    for(int s=i+j; s<orbitals; s++) {
                        int l = s - i - j - k;
                        m_listOfStates(counter,0) = i;
                        m_listOfStates(counter,1) = j;
                        m_listOfStates(counter,2) = k;
                        m_listOfStates(counter,3) = l;
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
}

double H0(const double x) {return 1;}
double H1(const double x) {return 2*x;}
double H2(const double x) {return 4*pow(x, 2) - 2;}
double H3(const double x) {return 8*pow(x, 3) - 12*x;}
double H4(const double x) {return 16*pow(x, 4) - 48*pow(x, 2) + 12;}
double H5(const double x) {return 32*pow(x, 5) - 160*pow(x, 3) + 120*x;}
double H6(const double x) {return 64*pow(x, 6) - 480*pow(x, 4) + 720*pow(x, 2) - 120;}
double H7(const double x) {return 128*pow(x, 7) - 1344*pow(x, 5) + 3360*pow(x, 3) - 1680*x;}
double H8(const double x) {return 256*pow(x, 8) - 3584*pow(x, 6) + 13440*pow(x, 4) - 13440*pow(x, 2) + 1680;}
double H9(const double x) {return 512*pow(x, 9) - 9216*pow(x, 7) + 48384*pow(x, 5) - 80640*pow(x, 3) + 30240*x;}
double H10(const double x) {return 1024*pow(x, 10) - 23040*pow(x, 8) + 161280*pow(x, 6) - 403200*pow(x, 4) + 302400*pow(x, 2) - 30240;}
double H11(const double x) {return 2048*pow(x, 11) - 56320*pow(x, 9) + 506880*pow(x, 7) - 1774080*pow(x, 5) + 2217600*pow(x, 3) - 665280*x;}
double H12(const double x) {return 4096*pow(x, 12) - 135168*pow(x, 10) + 1520640*pow(x, 8) - 7096320*pow(x, 6) + 13305600*pow(x, 4) - 7983360*pow(x, 2) + 665280;}
double H13(const double x) {return 8192*pow(x, 13) - 319488*pow(x, 11) + 4392960*pow(x, 9) - 26357760*pow(x, 7) + 69189120*pow(x, 5) - 69189120*pow(x, 3) + 17297280*x;}
double H14(const double x) {return 16384*pow(x, 14) - 745472*pow(x, 12) + 12300288*pow(x, 10) - 92252160*pow(x, 8) + 322882560*pow(x, 6) - 484323840*pow(x, 4) + 242161920*pow(x, 2) - 17297280;}
double H15(const double x) {return 32768*pow(x, 15) - 1720320*pow(x, 13) + 33546240*pow(x, 11) - 307507200*pow(x, 9) + 1383782400*pow(x, 7) - 2905943040*pow(x, 5) + 2421619200*pow(x, 3) - 518918400*x;}
double H16(const double x) {return 65536*pow(x, 16) - 3932160*pow(x, 14) + 89456640*pow(x, 12) - 984023040*pow(x, 10) + 5535129600*pow(x, 8) - 15498362880*pow(x, 6) + 19372953600*pow(x, 4) - 8302694400*pow(x, 2) + 518918400;}
double H17(const double x) {return 131072*pow(x, 17) - 8912896*pow(x, 15) + 233963520*pow(x, 13) - 3041525760*pow(x, 11) + 20910489600*pow(x, 9) - 75277762560*pow(x, 7) + 131736084480*pow(x, 5) - 94097203200*pow(x, 3) + 17643225600*x;}

double DH0(const double x) {return 0;}
double DH1(const double x) {return 2;}
double DH2(const double x) {return 8*x;}
double DH3(const double x) {return 24*pow(x, 2) - 12;}
double DH4(const double x) {return 64*pow(x, 3) - 96*x;}
double DH5(const double x) {return 160*pow(x, 4) - 480*pow(x, 2) + 120;}
double DH6(const double x) {return 384*pow(x, 5) - 1920*pow(x, 3) + 1440*x;}
double DH7(const double x) {return 896*pow(x, 6) - 6720*pow(x, 4) + 10080*pow(x, 2) - 1680;}
double DH8(const double x) {return 2048*pow(x, 7) - 17920*pow(x, 5) + 53760*pow(x, 3) - 26880*x;}
double DH9(const double x) {return 4608*pow(x, 8) - 64512*pow(x, 6) + 241920*pow(x, 4) - 241920*pow(x, 2) + 30240;}
double DH10(const double x) {return 10240*pow(x, 9) - 184320*pow(x, 7) + 967680*pow(x, 5) - 1612800*pow(x, 3) + 604800*x;}
double DH11(const double x) {return 22528*pow(x, 10) - 506880*pow(x, 8) + 3548160*pow(x, 6) - 8870400*pow(x, 4) + 6652800*pow(x, 2) - 665280;}
double DH12(const double x) {return 49152*pow(x, 11) - 1351680*pow(x, 9) + 12165120*pow(x, 7) - 42577920*pow(x, 5) + 53222400*pow(x, 3) - 15966720*x;}
double DH13(const double x) {return 106496*pow(x, 12) - 3514368*pow(x, 10) + 39536640*pow(x, 8) - 184504320*pow(x, 6) + 345945600*pow(x, 4) - 207567360*pow(x, 2) + 17297280;}
double DH14(const double x) {return 229376*pow(x, 13) - 8945664*pow(x, 11) + 123002880*pow(x, 9) - 738017280*pow(x, 7) + 1937295360*pow(x, 5) - 1937295360*pow(x, 3) + 484323840*x;}
double DH15(const double x) {return 491520*pow(x, 14) - 22364160*pow(x, 12) + 369008640*pow(x, 10) - 2767564800*pow(x, 8) + 9686476800*pow(x, 6) - 14529715200*pow(x, 4) + 7264857600*pow(x, 2) - 518918400;}
double DH16(const double x) {return 1048576*pow(x, 15) - 55050240*pow(x, 13) + 1073479680*pow(x, 11) - 9840230400*pow(x, 9) + 44281036800*pow(x, 7) - 92990177280*pow(x, 5) + 77491814400*pow(x, 3) - 16605388800*x;}
double DH17(const double x) {return 2228224*pow(x, 16) - 133693440*pow(x, 14) + 3041525760*pow(x, 12) - 33456783360*pow(x, 10) + 188194406400*pow(x, 8) - 52694337920*pow(x, 6) + 658680422400*pow(x, 4) - 282291609600*pow(x, 2) + 17643225600;}

double Hermite::evaluate(const double x, const int n) {
    //Hermite polynomial of n'th degree
    bool hardcoded = true;
    if(hardcoded) {
        switch(n) {
            case 0: return H0(x);
            case 1: return H1(x);
            case 2: return H2(x);
            case 3: return H3(x);
            case 4: return H4(x);
            case 5: return H5(x);
            case 6: return H6(x);
            case 7: return H7(x);
            case 8: return H8(x);
            case 9: return H9(x);
            case 10: return H10(x);
            case 11: return H11(x);
            case 12: return H12(x);
            case 13: return H13(x);
            case 14: return H14(x);
            case 15: return H15(x);
            case 16: return H16(x);
            case 17: return H17(x);
            default: return 0;
            if (n > 17) {
                return -1;
            }
        }
    }
    else {
        if(n == 0) {
            return 1;
        }
        else if(n == 1) {
            return 2 * m_omegaSqrt * x;
        }
        else {
            return 2 * (m_omegaSqrt * x * evaluate(x,n-1) - (n-1) * evaluate(x,n-2));
        }
    }
}

double Hermite::evaluateDerivative(const double x, const int n) {
    //First derivative of Hermite polynomial of n'th degree
    bool hardcoded = true;

    if(hardcoded || n < 18) {
        switch(n) {
            case 0: return DH0(x);
            case 1: return DH1(x);
            case 2: return DH2(x);
            case 3: return DH3(x);
            case 4: return DH4(x);
            case 5: return DH5(x);
            case 6: return DH6(x);
            case 7: return DH7(x);
            case 8: return DH8(x);
            case 9: return DH9(x);
            case 10: return DH10(x);
            case 11: return DH11(x);
            case 12: return DH12(x);
            case 13: return DH13(x);
            case 14: return DH14(x);
            case 15: return DH15(x);
            case 16: return DH16(x);
            case 17: return DH17(x);
            default: return 0;
        }
    }
    else {
        if(n == 0) {
            return 0;
        }
        else {
            return 2 * m_omegaSqrt * n * evaluate(x,n-1);
        }
    }
}

double Hermite::evaluateSecondDerivative(const double x, const int n) {
    //Second derivative of Hermite polynomial of n'th degree
    if(n < 2) {
        return 0;
    }
    else {
        return 4 * m_omegaSqrt * n * (n-1) * evaluate(x,n-2);
    }
}

double Hermite::basisElement(const unsigned int n, const Eigen::VectorXd positions) {
    double prod = 1;
    for(unsigned short i=0; i<m_numberOfDimensions; i++) {
        prod *= evaluate(positions(i), m_listOfStates(n, i));
    }
    return prod;
}

double Hermite::basisElementDer(const unsigned int n, const unsigned int i, const Eigen::VectorXd positions) {
    // i is the dimension we are derivating with respect to
    double prod = evaluateDerivative(positions(i), m_listOfStates(n, i));
    for(unsigned short j=0; j<m_numberOfDimensions; j++) {
        if(i != j) {
            prod *= evaluate(positions(j), m_listOfStates(n, j));
        }
    }
    return prod;
}

double Hermite::basisElementSecDer(const unsigned int n, const unsigned int i, const Eigen::VectorXd positions) {
    // i is the dimension we are derivating with respect to
    double prod = evaluateSecondDerivative(positions(i), m_listOfStates(n, i));
    for(unsigned int j=0; j<m_numberOfDimensions; j++) {
        if(i != j) {
            prod *= evaluate(positions(j), m_listOfStates(n, j));
        }
    }
    return prod;
}



