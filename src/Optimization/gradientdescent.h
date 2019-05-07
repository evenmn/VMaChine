#pragma once
#include "optimization.h"
#include <Eigen/Dense>

class GradientDescent : public Optimization {
public:
    GradientDescent(System* system, const double gamma, const double monotonicExp);
    unsigned int    getNumberOfBatches() { return m_numberOfBatches; }
    Eigen::MatrixXd updateParameters();

private:
    unsigned int    m_iter              = 0;
    unsigned int    m_numberOfBatches   = 1;
    double          m_gamma             = 0;
    double          m_monotonicExp      = 0;
    Eigen::MatrixXd m_v;
};
