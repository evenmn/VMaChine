#pragma once
#include "basis.h"

class HydrogenOrbital : public Basis
{
public:
    HydrogenOrbital(System *system);

    double evaluate(double x, int n);
    double evaluateDerivative(double x, int n);
    double evaluateSecondDerivative(double x, int n);

    void initialize();
    void setParameters(Eigen::VectorXd parameters);

    double basisElement(const int n, Eigen::VectorXd position);
    double basisElementDer(const int n, const int i, Eigen::VectorXd position);
    double basisElementSecDer(const int n, const int i, Eigen::VectorXd position);
    double basisElementPar(const int n, Eigen::VectorXd position);

    double evaluateCart(Eigen::VectorXd position, int n, int l, int m);
    double evaluateCartDerivative(Eigen::VectorXd position, int i, int n, int l, int m);
    double evaluateCartSecondDerivative(Eigen::VectorXd position, int i, int n, int l, int m);

    std::string getLabel() { return m_label; }

    void generateLOS();
    void numberOfOrbitalss();

private:
    int m_Z = 1;
    int m_numberOfOrbitalss = 1;
    int m_numberOfShells = 1;
    double m_alpha = 1;

    std::string m_label = "hydrogen orbital";

    Eigen::MatrixX3i m_LOS;
};
