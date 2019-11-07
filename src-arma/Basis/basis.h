#pragma once
#include <armadillo>
#include <fstream>
#include <iostream>
#include <mpi.h>

class Basis
{
public:
    Basis(class System *system);
    void numberOfOrbitals();
    void generateListOfStates();

    virtual void initialize() = 0;
    virtual void setParameters(arma::vec parameters) = 0;

    virtual double evaluate(double x, arma::uword n) = 0;
    virtual double evaluateDerivative(double x, arma::uword n) = 0;
    virtual double evaluateSecondDerivative(double x, arma::uword n) = 0;

    virtual double basisElement(const arma::uword n, arma::vec positions) = 0;
    virtual double basisElementDer(const arma::uword n, const arma::uword i, arma::vec positions) = 0;
    virtual double basisElementSecDer(const arma::uword n, const arma::uword i, arma::vec positions) = 0;
    virtual double basisElementPar(const arma::uword n, arma::vec position) = 0;

    virtual ~Basis() = 0;

    unsigned long long factorial(const arma::uword n);
    arma::uword factorialDifference(const arma::uword high, const arma::uword low);
    arma::uword binomial(const arma::uword n, const arma::uword p);
    void writeFileContent(std::string fileName, arma::mat &matrix);
    arma::uword fileLength(std::string fileName);

protected:
    arma::uword m_numberOfParticles = 0;
    arma::uword m_numberOfDimensions = 0;
    arma::uword m_numberOfOrbitals = 0;

    arma::umat m_listOfStates;

    class System *m_system = nullptr;
};
