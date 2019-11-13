#pragma once
#include "../Eigen/Dense"
#include <iostream>

class Layer
{
public:
    Layer(class System *system);
    virtual void initialize() = 0;
    virtual Eigen::VectorXd evaluate(Eigen::VectorXd a0) = 0;
    virtual Eigen::VectorXd activate(Eigen::VectorXd a0) = 0;
    virtual Eigen::VectorXd activateDer(Eigen::VectorXd z) = 0;
    virtual Eigen::VectorXd activateSecDer(Eigen::VectorXd z) = 0;
    virtual Eigen::VectorXd calculateDelta(Eigen::VectorXd z) = 0;
    virtual Eigen::MatrixXd calculateGradient(Eigen::VectorXd z, Eigen::VectorXd a0) = 0;
    virtual Eigen::MatrixXd updateWeights(Eigen::VectorXd z, Eigen::VectorXd a0) = 0;

    virtual ~Layer() = 0;

protected:
    int m_numberOfElements = 0;
    int m_maxParameters = 0;

    int m_h0 = 0;
    int m_h1 = 0;

    Eigen::MatrixXd m_W;

    class System *m_system = nullptr;
    class Activation *m_activation = nullptr;
};
