#pragma once
#include "basis.h"

class None : public Basis {
public:
    None(System* system);
    //void numberOfOrbitals();
    //void generateListOfStates(int orbitals);

    void setParameters(Eigen::VectorXd parameters);

    double evaluate(double x, int n);
    double evaluateDerivative(double x, int n);
    double evaluateSecondDerivative(double x, int n);

    double basisElement(const int n, Eigen::VectorXd positions);
    double basisElementDer(const int n, const int i, Eigen::VectorXd positions);
    double basisElementSecDer(const int n, const int i, Eigen::VectorXd positions);
    double basisElementPar(const int n, Eigen::MatrixXd positionBlock);
};
