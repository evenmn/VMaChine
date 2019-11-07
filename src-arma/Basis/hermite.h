#pragma once
#include "basis.h"

class Hermite : public Basis
{
public:
    Hermite(System *system);
    //void numberOfOrbitals();
    //void generateListOfStates(arma::uword orbitals);

    void initialize();
    void setParameters(arma::vec parameters);

    double evaluate(double x, arma::uword n);
    double evaluateDerivative(double x, arma::uword n);
    double evaluateSecondDerivative(double x, arma::uword n);

    double basisElement(const arma::uword n, arma::vec positions);
    double basisElementDer(const arma::uword n, const arma::uword i, arma::vec positions);
    double basisElementSecDer(const arma::uword n, const arma::uword i, arma::vec positions);
    double basisElementPar(const arma::uword n, arma::vec position);

private:
    double m_omega = 1;
    double m_omegaSqrt = 1;
    //arma::mat m_listOfStates;
};
