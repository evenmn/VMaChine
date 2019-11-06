#include "none.h"
#include "../system.h"
#include <iostream>

None::None(System *system)
    : Basis(system)
{
    m_system = system;
    numberOfOrbitals();
}
/*
void None::numberOfOrbitals() {
    m_numberOfOrbitals = 1;
}
*/

void None::setParameters(arma::vec parameters) {}

double None::evaluate(double x, arma::uword n)
{
    return 1;
}

double None::evaluateDerivative(double x, arma::uword n)
{
    return 0;
}

double None::evaluateSecondDerivative(double x, arma::uword n)
{
    return 0;
}

double None::basisElement(const arma::uword n, arma::vec positions)
{
    return 1;
}

double None::basisElementDer(const arma::uword n, const arma::uword i, arma::vec positions)
{
    // i is the dimension we are derivating with respect to
    return 0;
}

double None::basisElementSecDer(const arma::uword n, const arma::uword i, arma::vec positions)
{
    // i is the dimension we are derivating with respect to
    return 0;
}

double None::basisElementPar(const arma::uword n, arma::vec position)
{
    return 0;
}
/*
void None::generateListOfStates(arma::uword orbitals) {
}
*/
