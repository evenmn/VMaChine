#pragma once
#include "optimization.h"
#include <Eigen/Dense>

class SGD : public Optimization {
public:
    SGD(System* system, const double gamma, const double monotonicExp);
    int             getNumberOfBatches() { return m_numberOfBatches; }
    Eigen::MatrixXd getAllInstantGradients();
    Eigen::MatrixXd updateParameters();
    Eigen::MatrixXd getEnergyGradient();

private:
    int    m_step = 0;
    int    m_numberOfBatches = 10;
    double m_gamma = 0;
    double m_monotonicExp = 0;
    Eigen::MatrixXd m_v;
};
