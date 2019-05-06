#pragma once
#include "basis.h"

class Hermite : public Basis {
public:
    Hermite(System* system);
    void numberOfOrbitals();
    void generateListOfStates(int orbitals);

    double evaluate(double x, int n);
    double evaluateDerivative(double x, int n);
    double evaluateSecondDerivative(double x, int n);

    double basisElement(const int n, Eigen::VectorXd positions);
    double basisElementDer(const int n, const int i, Eigen::VectorXd positions);
    double basisElementSecDer(const int n, const int i, Eigen::VectorXd positions);

private:
    double m_omega = 1;
    double m_omegaSqrt = 1;
    Eigen::MatrixXd m_listOfStates;
};
