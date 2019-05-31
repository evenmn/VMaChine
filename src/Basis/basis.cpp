#include "basis.h"

Basis::Basis(System *system) {
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

Basis::~Basis() {}
