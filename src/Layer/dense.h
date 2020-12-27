#pragma once
#include "layer.h"

class Dense : public Layer
{
public:
    Dense(class System *system, int numberOfUnits, class Activation *activation);
    int getNumberOfUnits() { return m_numberOfUnits; }
    Eigen::VectorXd getA() { return m_a; }
    Eigen::VectorXd getDA() { return m_da; }
    Eigen::VectorXd getDDA() { return m_dda; }
    Eigen::MatrixXd getWeights() { return m_W; }
    class Activation *getActivation() { return m_activation; }
    void initialize(int numberOfUnitsInPreviousLayer);
    Eigen::VectorXd evaluate(Eigen::VectorXd a0);
    Eigen::VectorXd activate();
    Eigen::VectorXd activateDer();
    Eigen::VectorXd activateSecDer();
    Eigen::VectorXd calculateDelta(Eigen::VectorXd delta1);
    void updateWeights(Eigen::MatrixXd m_WNew);
    Vector2l getWeightDim();
    Eigen::MatrixXd calculateParameterGradient();

private:
    int m_numberOfUnits = 0;
    int m_numberOfUnitsInPreviousLayer = 0;

    Eigen::VectorXd m_z;
    Eigen::VectorXd m_a0;
    Eigen::VectorXd m_a;
    Eigen::VectorXd m_da;
    Eigen::VectorXd m_dda;
    Eigen::VectorXd m_delta;

    Eigen::MatrixXd m_W;

    class Activation *m_activation = nullptr;
};
