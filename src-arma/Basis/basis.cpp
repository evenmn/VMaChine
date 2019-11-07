#include "basis.h"

Basis::Basis(System *system)
{
    m_system = system;
}

unsigned long long Basis::factorial(const arma::uword n)
{
    return (n <= 1) ? 1 : factorial(n - 1) * n;
}

arma::uword Basis::factorialDifference(const arma::uword high, const arma::uword low)
{
    return (high <= 1 || high <= low) ? 1 : factorialDifference(high - 1, low) * high;
}

arma::uword Basis::binomial(const arma::uword n, const arma::uword p)
{
    //Binomial coefficients, equal to magic numbers
    //return arma::uword(factorial(n+p)/(factorial(n)*factorial(p)));
    return factorialDifference(n + p, n) / factorial(p);
}

arma::uword Basis::fileLength(std::string fileName)
{
    std::ifstream inFile(fileName.c_str());
    std::ifstream::pos_type length =  std::count(std::istreambuf_iterator<char>(inFile),
                      std::istreambuf_iterator<char>(),
                      '\n');
    return arma::uword(length);
}

void Basis::writeFileContent(std::string fileName, arma::mat &matrix)
{
    std::ifstream inFile(fileName);
    if (inFile.is_open()) {
        std::string line;
        double value;
        arma::uword i = 0;
        while (std::getline(inFile, line)) {
            std::istringstream iss(line);
            arma::uword j = 0;
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
    arma::uword counter = 0;
    while (true) {
        arma::uword orb = 2 * Basis::binomial(counter, m_numberOfDimensions);
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
    arma::uword numberOfStates = Basis::binomial(m_numberOfOrbitals - 1, m_numberOfDimensions);
    m_listOfStates.zeros(numberOfStates, m_numberOfDimensions);
    arma::uword counter = 0;
    // Two dimensions
    if (m_numberOfDimensions == 2) {
        for (arma::uword i = 0; i < m_numberOfOrbitals; i++) {
            for (arma::uword j = 0; j < i + 1; j++) {
                m_listOfStates(counter, 0) = i - j;
                m_listOfStates(counter, 1) = j;
                counter += 1;
            }
        }
    }
    // Three dimensions
    else if (m_numberOfDimensions == 3) {
        for (arma::uword i = 0; i < m_numberOfOrbitals; i++) {
            for (arma::uword j = 0; j < i + 1; j++) {
                for (arma::uword k = 0; k < i - j + 1; k++) {
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
