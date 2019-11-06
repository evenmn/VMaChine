#pragma once
#include "basis.h"

class HermiteSpin : public Basis {
public:
    HermiteSpin(System* system);
    void numberOfOrbitals();
    void generateListOfStates(int orbitals);

    double evaluate(double x, int n);
    double evaluateDerivative(double x, int n);
    double evaluateSecondDerivative(double x, int n);

    double sphericalHarmonics(int l, int m, arma::vec positions);
    double sphericalHarmonicsDer(int l, int m, int i, arma::vec positions);

    double basisElement(const int n, arma::vec positions);
    double basisElementDer(const int n, const int i, arma::vec positions);
    double basisElementSecDer(const int n, const int i, arma::vec positions);

private:
    double m_omega = 1;
    double m_omegaSqrt = 1;
    arma::uvec m_listOfStates;
};
