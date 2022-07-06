#pragma once
#include "basis.h"

class Hermite : public Basis
{
public:
    Hermite(System *system);
    //void numberOfOrbitals();
    //void generateListOfStates(int orbitals);

    void initialize() override;
    void setParameters(Eigen::VectorXd parameters) override;

    double evaluate(double x, int n) override;
    double evaluateDerivative(double x, int n) override;
    double evaluateSecondDerivative(double x, int n) override;

    double basisElement(const int n, Eigen::VectorXd positions) override;
    double basisElementDer(const int n, const int i, Eigen::VectorXd positions) override;
    double basisElementSecDer(const int n, const int i, Eigen::VectorXd positions) override;
    double basisElementPar(const int n, Eigen::VectorXd position) override;

    std::string getLabel() override { return m_label; }

private:
    double m_omega = 1;
    double m_omegaSqrt = 1;
    //Eigen::MatrixXi m_listOfStates;
    std::string m_label = "Hermite";
};
