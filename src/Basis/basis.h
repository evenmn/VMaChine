#pragma once
#include <Eigen/Dense>
#include <mpi.h>

class Basis {
public:
    Basis(class System* system);
    virtual void numberOfOrbitals() = 0;
    virtual void generateListOfStates(int orbitals) = 0;

    virtual double basisElement(const int n, Eigen::VectorXd positions) = 0;
    virtual double basisElementDer(const int n, const int i, Eigen::VectorXd positions) = 0;
    virtual double basisElementSecDer(const int n, const int i, Eigen::VectorXd positions) = 0;

    virtual ~Basis() = 0;

    long long factorial(const int n);
    int factorialDifference(const int high, const int low);
    int binomial(const int n, const int p);

    int getNumberOfOrbitals() { return m_numberOfOrbitals; }

protected:
    class System* m_system = nullptr;
    int m_numberOfParticles = 0;
    int m_numberOfDimensions = 0;
    int m_numberOfOrbitals = 0;
};
