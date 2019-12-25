#pragma once
#include "layer.h"

class Input : public Layer
{
public:
    Input(class System *system, int numberOfUnits);
    int getNumberOfUnits() { return m_numberOfUnits; }
    Eigen::VectorXd getA() { return m_a; }
    Eigen::VectorXd getDA() { return m_da; }
    Eigen::VectorXd getDDA() { return m_dda; }
    Eigen::MatrixXd getWeights() { return m_W; }
    class Activation *getActivation() { return m_activation; }
    Eigen::VectorXd evaluate(Eigen::VectorXd a0);
    Eigen::VectorXd activate();
    Eigen::VectorXd activateDer();
    Eigen::VectorXd activateSecDer();
    Eigen::VectorXd calculateDelta(Eigen::VectorXd delta1);
    Eigen::MatrixXd calculateGradient();
    void updateWeights(Eigen::MatrixXd m_WNew);
    void initialize(int numberOfUnitsInPreviousLayer, double factor=0.01);
    Vector2l getWeightDim();

private:
    int m_numberOfUnits = 0;

    Eigen::VectorXd m_x;
    Eigen::VectorXd m_delta;
    Eigen::VectorXd m_a;
    Eigen::VectorXd m_da;
    Eigen::VectorXd m_dda;

    Eigen::MatrixXd m_W;

    class Activation *m_activation = nullptr;
};
