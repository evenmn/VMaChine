#pragma once
#include "basis.h"

class HartreeFockHermite : public Basis {
public:
    HartreeFockHermite(System* system);
    void numberOfOrbitals();
    double evaluate(double x, int n);
    double evaluateDerivative(double x, int n);
    double evaluateSecondDerivative(double x, int n);
    Eigen::MatrixXd generateListOfStates();

    std::string generateFileName();
    void readCoefficientFile();

private:
    double  m_omega = 1;
    double  m_omegaSqrt = 1;
    int     m_basisSize = 1;
    std::string m_path = "Path is not given yet";
    Eigen::MatrixXd m_coefficients;
    class Basis* m_hermite = nullptr;
};
