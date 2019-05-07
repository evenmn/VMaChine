#pragma once
#include <Eigen/Dense>

class Metropolis {
public:
    Metropolis(class System* system);
    virtual Eigen::VectorXd updatePositions()         { return m_positions; }
    virtual Eigen::VectorXd updateRadialVector()      { return m_radialVector; }
    virtual Eigen::MatrixXd updateDistanceMatrix()    { return m_distanceMatrix; }

    virtual bool            acceptMove() = 0;
    virtual                 ~Metropolis() = 0;

    double calculateDistanceMatrixElement   (const unsigned int i, const unsigned int j);
    void   calculateDistanceMatrixCross     (const unsigned int particle);
    double calculateRadialVectorElement     (const unsigned int particle);

protected:
    double          m_stepLength                = 0;
    double          m_diff                      = 0.5;

    unsigned short  m_numberOfDimensions        = 0;
    unsigned int    m_numberOfParticles         = 0;
    unsigned int    m_numberOfFreeDimensions    = 0;


    bool            m_calculateDistanceMatrix   = false;
    bool            m_calculateRadialVector     = false;

    Eigen::VectorXd m_positions;
    Eigen::VectorXd m_positionsOld;
    Eigen::VectorXd m_radialVector;
    Eigen::VectorXd m_radialVectorOld;
    Eigen::MatrixXd m_distanceMatrix;
    Eigen::MatrixXd m_distanceMatrixOld;

    class System*   m_system = nullptr;
};
