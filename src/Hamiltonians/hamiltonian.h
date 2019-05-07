#pragma once
#include <Eigen/Dense>
#include <mpi.h>

class Hamiltonian {
public:
    Hamiltonian(class System* system);
    virtual unsigned int    getGlobalArrayNeed() = 0;
    virtual double          computeLocalEnergy() = 0;
    virtual                 ~Hamiltonian()       = 0;

protected:
    unsigned int    m_numberOfParticles  = 0;
    unsigned short  m_numberOfDimensions = 0;

    bool            m_interaction        = true;

    Eigen::VectorXd m_positions;
    Eigen::VectorXd m_radialVector;
    Eigen::MatrixXd m_distanceMatrix;
    Eigen::MatrixXd m_parameters;

    class System* m_system = nullptr;
};

