#pragma once
#include <Eigen/Dense>

class Metropolis {
public:
    Metropolis(class System* system);
    Eigen::VectorXd updatePositions()      { return m_positions; }
    virtual bool acceptMove() = 0;
    virtual ~Metropolis() = 0;

protected:
    class System* m_system = nullptr;
    Eigen::VectorXd m_positions;
    int m_numberOfFreeDimensions = 0;
    double m_stepLength = 0;
    double m_diff = 0.5;
};
