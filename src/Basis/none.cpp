#include "none.h"
#include "../system.h"
#include <iostream>

None::None(System *system)  :
    Basis(system) {
    m_system                = system;
    numberOfOrbitals();
}

void None::numberOfOrbitals() {
    m_numberOfOrbitals = 1;
}

double None::basisElement(const unsigned int n, const Eigen::VectorXd positions) {
    return 1;
}

double None::basisElementDer(const unsigned int n, const unsigned int i, const Eigen::VectorXd positions) {
    // i is the dimension we are derivating with respect to
    return 0;
}

double None::basisElementSecDer(const unsigned int n, const unsigned int i, const Eigen::VectorXd positions) {
    // i is the dimension we are derivating with respect to
    return 0;
}

void None::generateListOfStates(const int orbitals) {
}
