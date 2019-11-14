#pragma once
#include "layer.h"

class Dense : public Layer
{
public:
    Dense(class System *system, int h_0, int h_1, class Activation *activation);
    void initialize();
    Eigen::VectorXd evaluate(Eigen::VectorXd a0);
    Eigen::VectorXd activate(Eigen::VectorXd a0);
    Eigen::VectorXd activateDer(Eigen::VectorXd z);
    Eigen::VectorXd activateSecDer(Eigen::VectorXd z);
    Eigen::VectorXd calculateDelta(Eigen::VectorXd z);
    Eigen::MatrixXd calculateGradient(Eigen::VectorXd z, Eigen::VectorXd a0);
    void updateWeights(Eigen::MatrixXd m_WNew);
    Eigen::Vector2i getWeightDim();

private:
    typedef Eigen::Matrix<long, 2, 1> Vector2l;

    int m_numberOfElements = 0;
    int m_maxParameters = 0;

    int m_h0 = 0;
    int m_h1 = 0;

    Eigen::MatrixXd m_W;

    class Activation *m_activation = nullptr;
};
