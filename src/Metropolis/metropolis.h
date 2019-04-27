#pragma once
#include <Eigen/Dense>

class Metropolis {
public:
    Metropolis(class System* system);
    virtual Eigen::VectorXd updatePositions()         { return m_positions; }
    virtual Eigen::VectorXd updateRadialVector()      { return m_radialVector; }
    virtual Eigen::MatrixXd updateDistanceMatrix()    { return m_distanceMatrix; }

    virtual bool acceptMove() = 0;
    virtual double calculateDistanceMatrixElement(int i, int j) = 0;
    virtual void   calculateDistanceMatrixCross(int particle) = 0;
    virtual double calculateRadialVectorElement(int particle) = 0;
    virtual ~Metropolis() = 0;

protected:
    class System*   m_system = nullptr;
    Eigen::VectorXd m_positions;
    Eigen::VectorXd m_positionsOld;
    Eigen::VectorXd m_radialVector;
    Eigen::VectorXd m_radialVectorOld;
    Eigen::MatrixXd m_distanceMatrix;
    Eigen::MatrixXd m_distanceMatrixOld;
    int             m_numberOfFreeDimensions = 0;
    int             m_numberOfDimensions = 0;
    int             m_numberOfParticles = 0;
    double          m_stepLength = 0;
    double          m_diff = 0.5;

    bool            m_calculateDistanceMatrix = false;
    bool            m_calculateRadialVector = false;
};
