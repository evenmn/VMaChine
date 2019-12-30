#pragma once
#include <Eigen/Dense>
#include <iostream>

class Layer
{
public:
    typedef Eigen::Matrix<long, 2, 1> Vector2l;

    Layer(class System *system);
    virtual int getNumberOfUnits() = 0;
    virtual Eigen::VectorXd getA() = 0;
    virtual Eigen::VectorXd getDA() = 0;
    virtual Eigen::VectorXd getDDA() = 0;
    virtual Eigen::MatrixXd getWeights() = 0;
    virtual class Activation *getActivation() = 0;
    virtual Vector2l getWeightDim() = 0;
    virtual void initialize(int numberOfUnitsInPreviousLayer) = 0;
    virtual Eigen::VectorXd evaluate(Eigen::VectorXd a0) = 0;
    virtual Eigen::VectorXd activate() = 0;
    virtual Eigen::VectorXd activateDer() = 0;
    virtual Eigen::VectorXd activateSecDer() = 0;
    virtual Eigen::VectorXd calculateDelta(Eigen::VectorXd delta1) = 0;
    virtual void updateWeights(Eigen::MatrixXd WNew) = 0;

    virtual ~Layer() = 0;



protected:
    class System *m_system = nullptr;
};
