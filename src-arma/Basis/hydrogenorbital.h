#pragma once
#include "basis.h"

class HydrogenOrbital : public Basis
{
public:
    HydrogenOrbital(System *system);

    double evaluate(double x, arma::uword n);
    double evaluateDerivative(double x, arma::uword n);
    double evaluateSecondDerivative(double x, arma::uword n);

    void initialize();
    void setParameters(arma::vec parameters);

    double basisElement(const arma::uword n, arma::vec position);
    double basisElementDer(const arma::uword n, const arma::uword i, arma::vec position);
    double basisElementSecDer(const arma::uword n, const arma::uword i, arma::vec position);
    double basisElementPar(const arma::uword n, arma::vec position);

    double evaluateCart(arma::vec position, arma::uword n, arma::uword l, arma::sword m);
    double evaluateCartDerivative(arma::vec position, arma::uword i, arma::uword n, arma::uword l, arma::sword m);
    double evaluateCartSecondDerivative(arma::vec position, arma::uword i, arma::uword n, arma::uword l, arma::sword m);

    void generateLOS();
    void numberOfOrbitalss();

private:
    arma::uword m_Z = 1;
    arma::uword m_numberOfOrbitalss = 1;
    arma::uword m_numberOfShells = 1;
    double m_alpha = 1;

    arma::imat m_LOS;
};
