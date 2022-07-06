#pragma once
#include "basis.h"

class HartreeFock : public Basis
{
public:
    HartreeFock(System *system, Basis *basis);
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
    std::string generateFileName();
    void readCoefficientFile();

private:
    double m_omega = 1;
    double m_omegaSqrt = 1;
    long m_basisSize = 1;
    std::string m_path = "Path is not given yet";
    std::string m_label = "Hartree-Fock";

    Eigen::MatrixXd m_coefficients;
    //Eigen::MatrixXi m_listOfStates;
    class Basis *m_basis = nullptr;
};
