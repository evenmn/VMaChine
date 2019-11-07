#pragma once
#include "basis.h"

class HartreeFock : public Basis
{
public:
    HartreeFock(System *system, Basis *basis);
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

    std::string generateFileName();
    void readCoefficientFile();

private:
    double m_omega = 1;
    double m_omegaSqrt = 1;
    arma::uword m_basisSize = 1;
    std::string m_path = "Path is not given yet";
    arma::mat m_coefficients;
    //arma::mat m_listOfStates;
    class Basis *m_basis = nullptr;
};
