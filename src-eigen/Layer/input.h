#pragma once
#include "layer.h"

class Input : public Layer
{
public:
    Input(class System *system, int numberOfUnits);
    int getNumberOfUnits() { return m_numberOfUnits; }
    Eigen::VectorXd evaluate();
    Eigen::VectorXd activate(Eigen::VectorXd a0);
    Eigen::VectorXd calculateDelta(Eigen::VectorXd delta0);
    Eigen::MatrixXd calculateGradient();
    void updateWeights(Eigen::MatrixXd m_WNew);
    void initialize(int numberOfUnitsInPreviousLayer);
    Vector2l getWeightDim();

private:
    int m_numberOfUnits = 0;

    Eigen::VectorXd m_x;
    Eigen::VectorXd m_delta;
};
