#pragma once
#include "optimization.h"
#include <Eigen/Dense>

class GradientDescent : public Optimization {
public:
    GradientDescent(System* system, const double gamma, const double monotonicExp);
    int             getNumberOfBatches() { return m_numberOfBatches; }
    Eigen::MatrixXd updateParameters();
    Eigen::MatrixXd getEnergyGradient();

private:
    int    m_step = 0;
    int    m_numberOfBatches = 1;
    double m_gamma = 0;
    double m_monotonicExp = 0;
    Eigen::MatrixXd m_v;
};
