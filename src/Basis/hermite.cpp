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
        int orb = 2 * Basis::binomial(counter, m_numberOfDimensions);
        if(orb == m_numberOfParticles) {
            m_numberOfOrbitals = counter+1;
            break;
        }
        else if(orb > m_numberOfParticles) {
            std::cout << "This program supports closed-shells only. Please choose a number "
                         "of particles such that the orbital is full" << std::endl;
            MPI_Finalize();
            exit(0);
        }
        counter += 1;
    }
}

void Hermite::generateListOfStates(int orbitals) {
    // Returns the index list used in Slater
    // For instance (0,0), (1,0), (0,1) for 6P in 2D
    //              (0,0,0), (1,0,0), (0,1,0), (0,0,1) for 8P in 3D etc..
    int numberOfStates = Basis::binomial(orbitals-1, m_numberOfDimensions);
    m_listOfStates = Eigen::MatrixXi::Zero(numberOfStates, m_numberOfDimensions);
    int counter = 0;
    // Two dimensions
    if (m_numberOfDimensions == 2) {
        for(int i=0; i<orbitals; i++) {
            for(int j=0; j<i+1; j++) {
                m_listOfStates(counter,0) = i-j;
                m_listOfStates(counter,1) = j;
                counter += 1;
            }
        }
    }
    // Three dimensions
    else if (m_numberOfDimensions == 3) {
        for(int i=0; i<orbitals; i++) {
            for(int j=0; j<i+1; j++) {
                for(int k=0; k<i-j+1; k++) {
                    m_listOfStates(counter,0) = i-j-k;
                    m_listOfStates(counter,1) = j;
                    m_listOfStates(counter,2) = k;
                    counter += 1;
                }
            }
        }
    }
    else {
        std::cout << "Number of dimensions should be either 2 or 3" << std::endl;
        exit(0);
    }
}

double Hermite::basisElement(const int n, Eigen::VectorXd positions) {
    double prod = 1;
    for(int i=0; i<m_numberOfDimensions; i++) {
        prod *= evaluate(positions(i), int(m_listOfStates(n, i)));
    }
    return prod;
}

double Hermite::basisElementDer(const int n, const int i, Eigen::VectorXd positions) {
    // i is the dimension we are derivating with respect to
    double prod = evaluateDerivative(positions(i), m_listOfStates(n, i));
    for(int j=0; j<m_numberOfDimensions; j++) {
        if(i != j) {
            prod *= evaluate(positions(j), m_listOfStates(n, j));
        }
    }
    return prod;
}

double Hermite::basisElementSecDer(const int n, const int i, Eigen::VectorXd positions) {
    // i is the dimension we are derivating with respect to
    double prod = evaluateSecondDerivative(positions(i), m_listOfStates(n, i));
    for(int j=0; j<m_numberOfDimensions; j++) {
        if(i != j) {
            prod *= evaluate(positions(j), m_listOfStates(n, j));
        }
    }
    return prod;
}

double HH0(const double x) {return 1;}
double HH1(const double x) {return 2*x;}
double HH2(const double x) {return 4*pow(x, 2) - 2;}
double HH3(const double x) {return 8*pow(x, 3) - 12*x;}
double HH4(const double x) {return 16*pow(x, 4) - 48*pow(x, 2) + 12;}
double HH5(const double x) {return 32*pow(x, 5) - 160*pow(x, 3) + 120*x;}
double HH6(const double x) {return 64*pow(x, 6) - 480*pow(x, 4) + 720*pow(x, 2) - 120;}
double HH7(const double x) {return 128*pow(x, 7) - 1344*pow(x, 5) + 3360*pow(x, 3) - 1680*x;}
double HH8(const double x) {return 256*pow(x, 8) - 3584*pow(x, 6) + 13440*pow(x, 4) - 13440*pow(x, 2) + 1680;}
double HH9(const double x) {return 512*pow(x, 9) - 9216*pow(x, 7) + 48384*pow(x, 5) - 80640*pow(x, 3) + 30240*x;}
double HH10(const double x) {return 1024*pow(x, 10) - 23040*pow(x, 8) + 161280*pow(x, 6) - 403200*pow(x, 4) + 302400*pow(x, 2) - 30240;}
double HH11(const double x) {return 2048*pow(x, 11) - 56320*pow(x, 9) + 506880*pow(x, 7) - 1774080*pow(x, 5) + 2217600*pow(x, 3) - 665280*x;}
double HH12(const double x) {return 4096*pow(x, 12) - 135168*pow(x, 10) + 1520640*pow(x, 8) - 7096320*pow(x, 6) + 13305600*pow(x, 4) - 7983360*pow(x, 2) + 665280;}
double HH13(const double x) {return 8192*pow(x, 13) - 319488*pow(x, 11) + 4392960*pow(x, 9) - 26357760*pow(x, 7) + 69189120*pow(x, 5) - 69189120*pow(x, 3) + 17297280*x;}
double HH14(const double x) {return 16384*pow(x, 14) - 745472*pow(x, 12) + 12300288*pow(x, 10) - 92252160*pow(x, 8) + 322882560*pow(x, 6) - 484323840*pow(x, 4) + 242161920*pow(x, 2) - 17297280;}
double HH15(const double x) {return 32768*pow(x, 15) - 1720320*pow(x, 13) + 33546240*pow(x, 11) - 307507200*pow(x, 9) + 1383782400*pow(x, 7) - 2905943040*pow(x, 5) + 2421619200*pow(x, 3) - 518918400*x;}
double HH16(const double x) {return 65536*pow(x, 16) - 3932160*pow(x, 14) + 89456640*pow(x, 12) - 984023040*pow(x, 10) + 5535129600*pow(x, 8) - 15498362880*pow(x, 6) + 19372953600*pow(x, 4) - 8302694400*pow(x, 2) + 518918400;}
double HH17(const double x) {return 131072*pow(x, 17) - 8912896*pow(x, 15) + 233963520*pow(x, 13) - 3041525760*pow(x, 11) + 20910489600*pow(x, 9) - 75277762560*pow(x, 7) + 131736084480*pow(x, 5) - 94097203200*pow(x, 3) + 17643225600*x;}

double DHH0(const double x) {return 0;}
double DHH1(const double x) {return 2;}
double DHH2(const double x) {return 8*x;}
double DHH3(const double x) {return 24*pow(x, 2) - 12;}
double DHH4(const double x) {return 64*pow(x, 3) - 96*x;}
double DHH5(const double x) {return 160*pow(x, 4) - 480*pow(x, 2) + 120;}
double DHH6(const double x) {return 384*pow(x, 5) - 1920*pow(x, 3) + 1440*x;}
double DHH7(const double x) {return 896*pow(x, 6) - 6720*pow(x, 4) + 10080*pow(x, 2) - 1680;}
double DHH8(const double x) {return 2048*pow(x, 7) - 17920*pow(x, 5) + 53760*pow(x, 3) - 26880*x;}
double DHH9(const double x) {return 4608*pow(x, 8) - 64512*pow(x, 6) + 241920*pow(x, 4) - 241920*pow(x, 2) + 30240;}
double DHH10(const double x) {return 10240*pow(x, 9) - 184320*pow(x, 7) + 967680*pow(x, 5) - 1612800*pow(x, 3) + 604800*x;}
double DHH11(const double x) {return 22528*pow(x, 10) - 506880*pow(x, 8) + 3548160*pow(x, 6) - 8870400*pow(x, 4) + 6652800*pow(x, 2) - 665280;}
double DHH12(const double x) {return 49152*pow(x, 11) - 1351680*pow(x, 9) + 12165120*pow(x, 7) - 42577920*pow(x, 5) + 53222400*pow(x, 3) - 15966720*x;}
double DHH13(const double x) {return 106496*pow(x, 12) - 3514368*pow(x, 10) + 39536640*pow(x, 8) - 184504320*pow(x, 6) + 345945600*pow(x, 4) - 207567360*pow(x, 2) + 17297280;}
double DHH14(const double x) {return 229376*pow(x, 13) - 8945664*pow(x, 11) + 123002880*pow(x, 9) - 738017280*pow(x, 7) + 1937295360*pow(x, 5) - 1937295360*pow(x, 3) + 484323840*x;}
double DHH15(const double x) {return 491520*pow(x, 14) - 22364160*pow(x, 12) + 369008640*pow(x, 10) - 2767564800*pow(x, 8) + 9686476800*pow(x, 6) - 14529715200*pow(x, 4) + 7264857600*pow(x, 2) - 518918400;}
double DHH16(const double x) {return 1048576*pow(x, 15) - 55050240*pow(x, 13) + 1073479680*pow(x, 11) - 9840230400*pow(x, 9) + 44281036800*pow(x, 7) - 92990177280*pow(x, 5) + 77491814400*pow(x, 3) - 16605388800*x;}
double DHH17(const double x) {return 2228224*pow(x, 16) - 133693440*pow(x, 14) + 3041525760*pow(x, 12) - 33456783360*pow(x, 10) + 188194406400*pow(x, 8) - 52694337920*pow(x, 6) + 658680422400*pow(x, 4) - 282291609600*pow(x, 2) + 17643225600;}

double DDHH0(const double x) {return 0;}
double DDHH1(const double x) {return 0;}
double DDHH2(const double x) {return 8;}
double DDHH3(const double x) {return 48*x;}
double DDHH4(const double x) {return 192*pow(x, 2) - 96;}
double DDHH5(const double x) {return 640*pow(x, 3) - 960*x;}
double DDHH6(const double x) {return 1920*pow(x, 4) - 5760*pow(x, 2) + 1440;}
double DDHH7(const double x) {return 5376*pow(x, 5) - 26880*pow(x, 3) + 20160*x;}
double DDHH8(const double x) {return 14336*pow(x, 6) - 107520*pow(x, 4) + 161280*pow(x, 2) - 26880;}
double DDHH9(const double x) {return 36864*pow(x, 7) - 387072*pow(x, 5) + 967680*pow(x, 3) - 483840*x;}
double DDHH10(const double x) {return 92160*pow(x, 8) - 1290240*pow(x, 6) + 4838400*pow(x, 4) - 4838400*pow(x, 2) + 604800;}
double DDHH11(const double x) {return 225280*pow(x, 9) - 4055040*pow(x, 7) + 21288960*pow(x, 5) - 35481600*pow(x, 3) + 13305600*x;}
double DDHH12(const double x) {return 540672*pow(x, 10) - 12165120*pow(x, 8) + 85155840*pow(x, 6) - 212889600*pow(x, 4) + 159667200*pow(x, 2) - 15966720;}
double DDHH13(const double x) {return 1277952*pow(x, 11) - 35143680*pow(x, 9) + 316293120*pow(x, 7) - 1107025920*pow(x, 5) + 1383782400*pow(x, 3) - 415134720*x;}
double DDHH14(const double x) {return 2981888*pow(x, 12) - 98402304*pow(x, 10) + 1107025920*pow(x, 8) - 5166120960*pow(x, 6) + 9686476800*pow(x, 4) - 5811886080*pow(x, 2) + 484323840;}
double DDHH15(const double x) {return 6881280*pow(x, 13) - 268369920*pow(x, 11) + 3690086400*pow(x, 9) - 22140518400*pow(x, 7) + 58118860800*pow(x, 5) - 58118860800*pow(x, 3) + 14529715200*x;}
double DDHH16(const double x) {return 15728640*pow(x, 14) - 715653120*pow(x, 12) + 11808276480*pow(x, 10) - 88562073600*pow(x, 8) + 309967257600*pow(x, 6) - 464950886400*pow(x, 4) + 232475443200*pow(x, 2) - 16605388800;}
double DDHH17(const double x) {return 35651584*pow(x, 15) - 1871708160*pow(x, 13) + 36498309120*pow(x, 11) - 334567833600*pow(x, 9) + 1505555251200*pow(x, 7) - 316166027520*pow(x, 5) + 2634721689600*pow(x, 3) - 565583219200*x;}

double Hermite::evaluate(double x, int n) {
    //Hermite polynomial of n'th degree
    bool hardcoded = true;

    if(hardcoded) {
        switch(n) {
            case 0: return HH0(x);
            case 1: return HH1(x);
            case 2: return HH2(x);
            case 3: return HH3(x);
            case 4: return HH4(x);
            case 5: return HH5(x);
            case 6: return HH6(x);
            case 7: return HH7(x);
            case 8: return HH8(x);
            case 9: return HH9(x);
            case 10: return HH10(x);
            case 11: return HH11(x);
            case 12: return HH12(x);
            case 13: return HH13(x);
            case 14: return HH14(x);
            case 15: return HH15(x);
            case 16: return HH16(x);
            case 17: return HH17(x);
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

    //return std::hermite(unsigned(n),x);
}

double Hermite::evaluateDerivative(double x, int n) {
    //First derivative of Hermite polynomial of n'th degree
    bool hardcoded = true;

    if(hardcoded || n < 18) {
        switch(n) {
            case 0: return DHH0(x);
            case 1: return DHH1(x);
            case 2: return DHH2(x);
            case 3: return DHH3(x);
            case 4: return DHH4(x);
            case 5: return DHH5(x);
            case 6: return DHH6(x);
            case 7: return DHH7(x);
            case 8: return DHH8(x);
            case 9: return DHH9(x);
            case 10: return DHH10(x);
            case 11: return DHH11(x);
            case 12: return DHH12(x);
            case 13: return DHH13(x);
            case 14: return DHH14(x);
            case 15: return DHH15(x);
            case 16: return DHH16(x);
            case 17: return DHH17(x);
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
    bool hardcoded = true;

    if(hardcoded || n < 18) {
        switch(n) {
            case 0: return DDHH0(x);
            case 1: return DDHH1(x);
            case 2: return DDHH2(x);
            case 3: return DDHH3(x);
            case 4: return DDHH4(x);
            case 5: return DDHH5(x);
            case 6: return DDHH6(x);
            case 7: return DDHH7(x);
            case 8: return DDHH8(x);
            case 9: return DDHH9(x);
            case 10: return DDHH10(x);
            case 11: return DDHH11(x);
            case 12: return DDHH12(x);
            case 13: return DDHH13(x);
            case 14: return DDHH14(x);
            case 15: return DDHH15(x);
            case 16: return DDHH16(x);
            case 17: return DDHH17(x);
            default: return 0;
        }
    }
    else {
        if(n < 2) {
            return 0;
        }
        else {
            return 4 * m_omegaSqrt * n * (n-1) * evaluate(x,n-2);
        }
    }
}



