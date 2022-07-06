#pragma once
#include "basis.h"

class HartreeFockHermite : public Basis {
public:
    HartreeFockHermite(System* system, Basis *basis);
    void numberOfOrbitals();
    void generateListOfStates(int orbitals);

    double basisElement(const int n, Eigen::VectorXd positions);
    double basisElementDer(const int n, const int i, Eigen::VectorXd positions);
    double basisElementSecDer(const int n, const int i, Eigen::VectorXd positions);

    std::string getLabel() { return m_label; }
    std::string generateFileName();
    void readCoefficientFile();

private:
    double  m_omega = 1;
    double  m_omegaSqrt = 1;
    long    m_basisSize = 1;
    std::string m_path = "Path is not given yet";
    std::string m_label = "Hartree-Fock Hermite";
    Eigen::MatrixXd m_coefficients;

    class Basis* m_hermite = nullptr;
    class Basis* m_basis = nullptr;
};
