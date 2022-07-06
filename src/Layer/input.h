#pragma once
#include "layer.h"

class Input : public Layer
{
public:
    Input(class System *system, int numberOfUnits);
    int getNumberOfUnits() override { return m_numberOfUnits; }
    Eigen::VectorXd getA() override { return m_a; }
    Eigen::VectorXd getDA() override { return m_da; }
    Eigen::VectorXd getDDA() override { return m_dda; }
    Eigen::MatrixXd getWeights() override { return m_W; }
    class Activation *getActivation() override { return m_activation; }
    Eigen::VectorXd evaluate(Eigen::VectorXd x) override;
    Eigen::VectorXd activate() override;
    Eigen::VectorXd activateDer() override;
    Eigen::VectorXd activateSecDer() override;
    Eigen::VectorXd calculateDelta(Eigen::VectorXd delta1) override;
    void updateWeights(Eigen::MatrixXd m_WNew) override;
    void initialize(int numberOfUnitsInPreviousLayer) override;
    Vector2l getWeightDim() override;

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
