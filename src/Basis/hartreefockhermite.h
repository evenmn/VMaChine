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

private:
    double  m_omega = 1;
    double  m_omegaSqrt = 1;
    int     m_basisSize = 1;
    Eigen::MatrixXd m_coefficients;
    class Basis* m_hermite = nullptr;
};
