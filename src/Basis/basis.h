#pragma once
#include <Eigen/Dense>
#include <mpi.h>

class Basis {
public:
    Basis(class System* system);
    virtual void numberOfOrbitals() = 0;
    virtual double evaluate(double x, int n) = 0;
    virtual double evaluateDerivative(double x, int n) = 0;
    virtual double evaluateSecondDerivative(double x, int n) = 0;
    virtual Eigen::MatrixXd generateListOfStates() = 0;
    virtual ~Basis() = 0;

    int factorial(const int n);
    double binomial(const int n, const int p);

    int getNumberOfOrbitals() { return m_numberOfOrbitals; }

protected:
    class System* m_system = nullptr;
    int m_numberOfParticles = 0;
    int m_numberOfDimensions = 0;
    int m_numberOfOrbitals = 0;
};
