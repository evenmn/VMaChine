#pragma once
#include "optimization.h"
#include <Eigen/Dense>

class GradientDescent : public Optimization {
public:
    GradientDescent(System* system, const double gamma);
    Eigen::VectorXd getImmediateGradients(WaveFunction* waveFunction);
    Eigen::MatrixXd getAllImmediateGradients();
    Eigen::MatrixXd updateParameters();
    Eigen::MatrixXd getEnergyGradient();

private:
    double m_gamma = 0;
    int    m_step = 0;
    Eigen::MatrixXd m_v;
};
