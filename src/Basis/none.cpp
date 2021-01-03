#include "none.h"
#include "../system.h"
#include <iostream>

None::None(System *system)
    : Basis(system)
{
    /* Constructor of the None class.
     * This is used when simulating small atomic systems. */
    m_system = system;
    //numberOfOrbitals();
}
/*
void None::numberOfOrbitals() {
    m_numberOfOrbitals = 1;
}
*/

void None::initialize() {}

void None::setParameters(Eigen::VectorXd /*parameters*/) {}

double None::evaluate(double /*x*/, int /*n*/)
{
    /* Evaluate None. */
    return 1;
}

double None::evaluateDerivative(double /*x*/, int /*n*/)
{
    return 0;
}

double None::evaluateSecondDerivative(double /*x*/, int /*n*/)
{
    return 0;
}

double None::basisElement(const int /*n*/, Eigen::VectorXd /*positions*/)
{
    return 1;
}

double None::basisElementDer(const int /*n*/, const int /*i*/, Eigen::VectorXd /*positions*/)
{
    // i is the dimension we are derivating with respect to
    return 0;
}

double None::basisElementSecDer(const int /*n*/, const int /*i*/, Eigen::VectorXd /*positions*/)
{
    // i is the dimension we are derivating with respect to
    return 0;
}

double None::basisElementPar(const int /*n*/, Eigen::VectorXd /*position*/)
{
    return 0;
}
/*
void None::generateListOfStates(int orbitals) {
}
*/
