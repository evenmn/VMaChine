#pragma once
#include "../Eigen/Dense"

class Hamiltonian {
public:
    Hamiltonian(class System* system);
    virtual int    getGlobalArrayNeed() = 0;
    virtual int    getNumberOfSources() = 0;
    virtual double getExternalEnergy()  = 0;
    virtual double computeLocalEnergy() = 0;
    virtual ~Hamiltonian() = 0;

    double getInteractionEnergy();

protected:
    int             m_numberOfParticles = 0;
    int             m_numberOfDimensions = 0;

    bool            m_interaction = false;

    Eigen::VectorXd m_positions;
    Eigen::VectorXd m_radialVector;
    Eigen::MatrixXd m_distanceMatrix;
    Eigen::MatrixXd m_parameters;

    class System*   m_system = nullptr;
};

