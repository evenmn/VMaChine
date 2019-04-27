#pragma once
#include <Eigen/Dense>

class InitialState {
public:
    InitialState(class System* system);
    virtual void setupInitialState() = 0;
    Eigen::VectorXd getParticles() { return m_positions; }
    Eigen::MatrixXd getDistanceMatrix() { return m_distanceMatrix; }
    Eigen::VectorXd getRadialVector() { return m_radialVector; }

    virtual double  calculateDistanceMatrixElement  (int i, int j) = 0;
    virtual void    calculateDistanceMatrix         () = 0;
    virtual double  calculateRadialVectorElement    (int particle) = 0;
    virtual void    calculateRadialVector           () = 0;

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

