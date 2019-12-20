#pragma once
#include <Eigen/Dense>
#include <iostream>

class Layer
{
public:
    typedef Eigen::Matrix<long, 2, 1> Vector2l;

    Layer(class System *system);
    virtual int getNumberOfUnits() = 0;
    virtual Vector2l getWeightDim() = 0;
    virtual void initialize(int numberOfUnitsInPreviousLayer) = 0;
    virtual Eigen::VectorXd activate(Eigen::VectorXd a0) = 0;
    virtual Eigen::VectorXd calculateDelta(Eigen::VectorXd delta0) = 0;
    virtual Eigen::MatrixXd calculateGradient() = 0;
    virtual void updateWeights(Eigen::MatrixXd WNew) = 0;

    virtual ~Layer() = 0;



protected:
    Eigen::MatrixXd m_W;

    class System *m_system = nullptr;
    class Activation *m_activation = nullptr;
};
