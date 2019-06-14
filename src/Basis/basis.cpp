#include "basis.h"

Basis::Basis(System *system) {
    m_system = system;
}

long long Basis::factorial(const int n) {
    return (n <= 1) ? 1 : factorial(n - 1) * n;
}

int Basis::factorialDifference(const int high, const int low) {
    return (high <= 1 || high <= low) ? 1 : factorialDifference(high - 1, low) * high;
}

int Basis::binomial(const int n, const int p) {
    //Binomial coefficients, equal to magic numbers
    //return int(factorial(n+p)/(factorial(n)*factorial(p)));
    return factorialDifference(n+p, n) / factorial(p);
}

std::ifstream::pos_type Basis::fileLength(std::string fileName) {
    std::ifstream inFile(fileName.c_str());
    return std::count(std::istreambuf_iterator<char>(inFile),
                      std::istreambuf_iterator<char>(), '\n');
}

void Basis::writeFileContentIntoEigenMatrix(std::string fileName, Eigen::MatrixXd &matrix) {
    std::ifstream inFile(fileName);
    if(inFile.is_open()) {
        std::string line;
        double value;
        int i = 0;
        while(std::getline(inFile, line)) {
            std::istringstream iss(line);
            int j = 0;
            while (iss >> value) {
                matrix(i,j) = value;
                j++;
            }
            i++;
        }
    }
    else {
        std::cout << "File '" << fileName << "' was not found" << std::endl;
        MPI_Finalize();
        exit(0);
    }
}

Eigen::MatrixXi Basis::generateListOfStates(int numberOfSources) {
    if(numberOfSources == 1) {
        int orbitals = 0;
        int i = 0;
        while(true) {
            int orb = 2 * binomial(i, m_numberOfDimensions);
            if(orb == m_numberOfParticles) {
                orbitals = i+1;
                break;
            }
            else if(orb > m_numberOfParticles) {
                std::cout << "Warning: An open shell configuration is chosen" << std::endl;
                orbitals = i+1;
                break;
            }
            i++;
        }

        // Returns the index list used in Slater
        // For instance (0,0), (1,0), (0,1) for 6P in 2D
        //              (0,0,0), (1,0,0), (0,1,0), (0,0,1) for 8P in 3D etc..
        int numberOfStates = binomial(orbitals-1, m_numberOfDimensions);
        Eigen::VectorXi listOfStates = Eigen::MatrixXi::Zero(numberOfStates, m_numberOfDimensions);
        int counter = 0;
        // Two dimensions
        if (m_numberOfDimensions == 2) {
            for(int i=0; i<orbitals; i++) {
                for(int j=0; j<i+1; j++) {
                    listOfStates(counter,0) = i-j;
                    listOfStates(counter,1) = j;
                    counter += 1;
                }
            }
            return listOfStates;
        }
        // Three dimensions
        else {
            for(int i=0; i<orbitals; i++) {
                for(int j=0; j<i+1; j++) {
                    for(int k=0; k<i-j+1; k++) {
                        listOfStates(counter,0) = i-j-k;
                        listOfStates(counter,1) = j;
                        listOfStates(counter,2) = k;
                        counter += 1;
                    }
                }
            }
            return listOfStates;
        }
    }
    else {
        int orbitalsLHS = 0;
        int orbitalsRHS = 0;
        int orb1 = 0;
        int orb2 = 0;
        int i = 0;
        while(true) {
            if(i==0) {
                orb1 = 0;
            }
            else {
                orb1 = 2 * binomial(i-1, m_numberOfDimensions);
            }
            int j = 0;
            while(j<2) {
                if(i==0) {
                    orb2 = 0;
                }
                else {
                    orb2 = 2 * binomial(j-1, m_numberOfDimensions);
                }
                if(orb1+orb2 == m_numberOfParticles) {
                    orbitalsLHS = i;
                    orbitalsRHS = j;
                    goto endloop;
                }
                else if(orb1 + orb2 > m_numberOfParticles) {
                    std::cout << "This program supports closed-shells only. Please choose a "
                                 "number of particles such that the orbital is full" << std::endl;
                    MPI_Finalize();
                    exit(0);
                }
                j++;
            }
            i++;
        }
        endloop:
        // Returns the index list used in Slater
        // For instance (0,0), (1,0), (0,1) for 6P in 2D
        //              (0,0,0), (1,0,0), (0,1,0), (0,0,1) for 8P in 3D etc..
        int numberOfStatesLHS = binomial(orbitalsLHS-1, m_numberOfDimensions);
        int numberOfStatesRHS = binomial(orbitalsRHS-1, m_numberOfDimensions);
        Eigen::MatrixXi listOfStates = Eigen::MatrixXi::Zero(numberOfStatesLHS+numberOfStatesRHS, m_numberOfDimensions);
        int counter = 0;
        // Two dimensions
        if (m_numberOfDimensions == 2) {
            for(int i=0; i<orbitalsLHS; i++) {
                for(int j=0; j<i+1; j++) {
                    listOfStates(counter,0) = i-j;
                    listOfStates(counter,1) = j;
                    counter += 1;
                }
            }
            for(int i=0; i<orbitalsRHS; i++) {
                for(int j=0; j<i+1; j++) {
                    listOfStates(counter,0) = i-j;
                    listOfStates(counter,1) = j;
                    counter += 1;
                }
            }
            return listOfStates;
        }
        // Three dimensions
        else {
            //std::cout << "jj" << std::endl;
            for(int i=0; i<orbitalsLHS; i++) {
                for(int j=0; j<i+1; j++) {
                    for(int k=0; k<i-j+1; k++) {
                        listOfStates(counter,0) = i-j-k;
                        listOfStates(counter,1) = j;
                        listOfStates(counter,2) = k;
                        counter += 1;
                    }
                }
            }
            for(int i=0; i<orbitalsRHS; i++) {
                for(int j=0; j<i+1; j++) {
                    for(int k=0; k<i-j+1; k++) {
                        listOfStates(counter,0) = i-j-k;
                        listOfStates(counter,1) = j;
                        listOfStates(counter,2) = k;
                        counter += 1;
                    }
                }
            }
            return listOfStates;
        }
    }
}

Basis::~Basis() {}
