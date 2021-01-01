#pragma once
#include <cassert>
#include <Eigen/Dense>

class Interaction
{
public:
    Interaction(class System *system);
    virtual int getGlobalArrayNeed() = 0;
    virtual std::string getLabel() = 0;

    virtual void initialize() = 0;
    virtual double getInteractionEnergy() = 0;
    virtual ~Interaction() = 0;

protected:
    int m_numberOfParticles = 0;
    int m_numberOfDimensions = 0;
    bool m_screening = false;

    double m_screeningStrength = 100;
    double m_dsl = 100;

    Eigen::VectorXd m_positions;
    Eigen::VectorXd m_radialVector;
    Eigen::MatrixXd m_distanceMatrix;
    Eigen::MatrixXd m_parameters;

    class System *m_system = nullptr;
};
