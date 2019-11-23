#pragma once
#include "layer.h"

class Dense : public Layer
{
public:
    Dense(class System *system, int h_0, int h_1, class Activation *activation);
    void initialize();
    Eigen::VectorXd evaluate();
    Eigen::VectorXd activate(Eigen::VectorXd a0);
    Eigen::VectorXd calculateDelta(Eigen::VectorXd delta0);
    Eigen::MatrixXd calculateGradient();
    void updateWeights(Eigen::MatrixXd m_WNew);
    Vector2l getWeightDim();

private:
    int m_h0 = 0;
    int m_h1 = 0;

    Eigen::VectorXd m_z;
    Eigen::VectorXd m_a;
    Eigen::VectorXd m_a0;
    Eigen::VectorXd m_delta;
};
