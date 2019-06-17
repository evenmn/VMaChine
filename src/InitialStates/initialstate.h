#pragma once
#include "../Eigen/Dense"

class InitialState {
public:
    InitialState(class System* system);
    virtual void setupInitialState() = 0;
    Eigen::VectorXd getParticles()      { return m_positions; }
    Eigen::MatrixXd getDistanceMatrix() { return m_distanceMatrix; }
    Eigen::VectorXd getRadialVector()   { return m_radialVector; }

    double  calculateDistanceMatrixElement  (int i, int j);
    void    calculateDistanceMatrix         ();
    double  calculateRadialVectorElement    (int particle);
    void    calculateRadialVector           ();

    virtual ~InitialState() = 0;

protected:
    class System* m_system = nullptr;
    int m_numberOfFreeDimensions = 0;
    int m_numberOfDimensions     = 0;
    int m_numberOfParticles      = 0;

    Eigen::VectorXd m_positions;
    Eigen::MatrixXd m_distanceMatrix;
    Eigen::VectorXd m_radialVector;
};

