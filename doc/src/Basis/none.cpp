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

double None::evaluate(double x, int n) {
    return 1;
}

double None::evaluateDerivative(double x, int n) {
    return 0;
}

double None::evaluateSecondDerivative(double x, int n) {
    return 0;
}

Eigen::MatrixXd None::generateListOfStates() {
    return Eigen::MatrixXd::Ones(1,1);
}
