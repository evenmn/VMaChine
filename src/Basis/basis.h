#pragma once
#include <Eigen/Dense>
#include <mpi.h>

class Basis {
public:
    Basis(class System* system);
    virtual void        numberOfOrbitals    () = 0;
    virtual void        generateListOfStates(const int orbitals) = 0;

    virtual double      basisElement        (const unsigned int n, const Eigen::VectorXd positions) = 0;
    virtual double      basisElementDer     (const unsigned int n, const unsigned int i, Eigen::VectorXd positions) = 0;
    virtual double      basisElementSecDer  (const unsigned int n, const unsigned int i, Eigen::VectorXd positions) = 0;

    virtual             ~Basis() = 0;

    unsigned long long  factorial           (const unsigned int n);
    unsigned int        binomial            (const unsigned int n, const unsigned int p);

    int                 getNumberOfOrbitals () { return m_numberOfOrbitals; }

protected:
    unsigned int        m_numberOfParticles  = 0;
    unsigned short      m_numberOfDimensions = 0;
    int                 m_numberOfOrbitals   = 0;

    class System*       m_system = nullptr;
};
