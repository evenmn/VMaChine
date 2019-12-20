#pragma once
#include "layer.h"

class Dense : public Layer
{
public:
    Dense(class System *system, int numberOfUnits, class Activation *activation);
    int getNumberOfUnits() { return m_numberOfUnits; }
    void initialize(int numberOfUnitsInPreviousLayer);
    Eigen::VectorXd evaluate();
    Eigen::VectorXd activate(Eigen::VectorXd a0);
    Eigen::VectorXd calculateDelta(Eigen::VectorXd delta0);
    Eigen::MatrixXd calculateGradient();
    void updateWeights(Eigen::MatrixXd m_WNew);
    Vector2l getWeightDim();

private:
    int m_numberOfUnits = 0;
    int m_numberOfUnitsInPreviousLayer = 0;

    Eigen::VectorXd m_z;
    Eigen::VectorXd m_a;
    Eigen::VectorXd m_a0;
    Eigen::VectorXd m_delta;


};
