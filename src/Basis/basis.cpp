#include "basis.h"

Basis::Basis(System *system) {}

long long Basis::factorial(const int n)
{
    return (n <= 1) ? 1 : factorial(n - 1) * n;
}

int Basis::factorialDifference(const int high, const int low)
{
    return (high <= 1 || high <= low) ? 1 : factorialDifference(high - 1, low) * high;
}

int Basis::binomial(const int n, const int p)
{
    //Binomial coefficients, equal to magic numbers
    //return int(factorial(n+p)/(factorial(n)*factorial(p)));
    return factorialDifference(n + p, n) / factorial(p);
}

std::ifstream::pos_type Basis::fileLength(std::string fileName)
{
    std::ifstream inFile(fileName.c_str());
    return std::count(std::istreambuf_iterator<char>(inFile),
                      std::istreambuf_iterator<char>(),
                      '\n');
}

void Basis::writeFileContentIntoEigenMatrix(std::string fileName, Eigen::MatrixXd &matrix)
{
    std::ifstream inFile(fileName);
    if (inFile.is_open()) {
        std::string line;
        double value;
        int i = 0;
        while (std::getline(inFile, line)) {
            std::istringstream iss(line);
            int j = 0;
            while (iss >> value) {
                matrix(i, j) = value;
                j++;
            }
            i++;
        }
    } else {
        std::cout << "File '" << fileName << "' was not found" << std::endl;
        MPI_Finalize();
        exit(0);
    }
}

void Basis::numberOfOrbitals()
{
    //Number of closed-shell orbitals
    int counter = 0;
    while (true) {
        int orb = 2 * Basis::binomial(counter, m_numberOfDimensions);
        if (orb == m_numberOfParticles) {
            m_numberOfOrbitals = counter + 1;
            break;
        } else if (orb > m_numberOfParticles) {
            std::cout << "Warning: An open shell is chosen" << std::endl;
            m_numberOfOrbitals = counter + 1;
            break;
        }
        counter += 1;
    }
}

void Basis::generateListOfStates()
{
    // Returns the index list used in Slater
    // For instance (0,0), (1,0), (0,1) for 6P in 2D
    //              (0,0,0), (1,0,0), (0,1,0), (0,0,1) for 8P in 3D etc..
    numberOfOrbitals();
    int numberOfStates = Basis::binomial(m_numberOfOrbitals - 1, m_numberOfDimensions);
    m_listOfStates = Eigen::MatrixXi::Zero(numberOfStates, m_numberOfDimensions);
    int counter = 0;
    // Two dimensions
    if (m_numberOfDimensions == 2) {
        for (int i = 0; i < m_numberOfOrbitals; i++) {
            for (int j = 0; j < i + 1; j++) {
                m_listOfStates(counter, 0) = i - j;
                m_listOfStates(counter, 1) = j;
                counter += 1;
            }
        }
    }
    // Three dimensions
    else if (m_numberOfDimensions == 3) {
        for (int i = 0; i < m_numberOfOrbitals; i++) {
            for (int j = 0; j < i + 1; j++) {
                for (int k = 0; k < i - j + 1; k++) {
                    m_listOfStates(counter, 0) = i - j - k;
                    m_listOfStates(counter, 1) = j;
                    m_listOfStates(counter, 2) = k;
                    counter += 1;
                }
            }
        }
    } else {
        std::cout << "Number of dimensions should be either 2 or 3" << std::endl;
        exit(0);
    }
}

Basis::~Basis() {}
