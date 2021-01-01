#pragma once
#include <Eigen/Dense>

class Hamiltonian
{
public:
    Hamiltonian(class System *system);
    virtual int getGlobalArrayNeed() = 0;
    virtual std::string getLabel() = 0;

    virtual void initialize() = 0;
    virtual double getExternalEnergy() = 0;
    virtual ~Hamiltonian() = 0;

protected:
    int m_numberOfParticles = 0;
    int m_numberOfDimensions = 0;

    Eigen::VectorXd m_positions;
    Eigen::VectorXd m_radialVector;
    Eigen::MatrixXd m_distanceMatrix;
    Eigen::MatrixXd m_parameters;

    class System *m_system = nullptr;
};
