#pragma once
#include "basis.h"

class HermiteExpansion : public Basis {
public:
    HermiteExpansion(System* system);
    void numberOfOrbitals();
    void generateListOfStates(int orbitals);

    double evaluate(double x, int n);
    double evaluateDerivative(double x, int n);
    double evaluateSecondDerivative(double x, int n);

    double basisElement(const int n, Eigen::VectorXd positions);
    double basisElementDer(const int n, const int i, Eigen::VectorXd positions);
    double basisElementSecDer(const int n, const int i, Eigen::VectorXd positions);

    std::string generateFileName();
    void readCoefficientFile();

private:
    double  m_omega = 1;
    double  m_omegaSqrt = 1;
    long    m_basisSize = 1;
    int     m_dim = 0;
    std::string m_path = "Path is not given yet";
    Eigen::MatrixXd m_coefficients;
    Eigen::MatrixXi m_listOfStates;
    class Basis* m_basis = nullptr;
};
