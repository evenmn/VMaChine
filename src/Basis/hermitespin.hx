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

    double sphericalHarmonics(int l, int m, Eigen::VectorXd positions);
    double sphericalHarmonicsDer(int l, int m, int i, Eigen::VectorXd positions);

    double basisElement(const int n, Eigen::VectorXd positions);
    double basisElementDer(const int n, const int i, Eigen::VectorXd positions);
    double basisElementSecDer(const int n, const int i, Eigen::VectorXd positions);

    std::string getLabel() { return m_label; }

private:
    double m_omega = 1;
    double m_omegaSqrt = 1;
    Eigen::MatrixXi m_listOfStates;
    std::string m_label = "Hermite spin";
};
