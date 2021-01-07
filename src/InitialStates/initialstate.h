#pragma once
#include <Eigen/Dense>
#include <cassert>
#include <iostream>

class InitialState
{
public:
    InitialState(class System *system);
    virtual std::string getLabel() = 0;
    virtual void setupInitialState() = 0;
    Eigen::VectorXd getParticles() { return m_positions; }
    Eigen::MatrixXd getDistanceMatrix() { return m_distanceMatrix; }
    Eigen::VectorXd getRadialVector() { return m_radialVector; }

    double calculateDistanceMatrixElement(int i, int j);
    void calculateDistanceMatrix();
    double calculateRadialVectorElement(int particle);
    void calculateRadialVector();

    virtual ~InitialState() = 0;

protected:
    class System *m_system = nullptr;
    int m_degreesOfFreedom = 0;
    int m_numberOfDimensions = 0;
    int m_numberOfParticles = 0;

    Eigen::VectorXd m_positions;
    Eigen::MatrixXd m_distanceMatrix;
    Eigen::VectorXd m_radialVector;
};
