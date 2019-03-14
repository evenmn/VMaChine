#pragma once
#include "optimization.h"
#include <Eigen/Dense>

class ASGD : public Optimization {
public:
    ASGD(System* system, const double gamma);
    int             getNumberOfBatches() { return m_numberOfBatches; }
    Eigen::VectorXd getInstantGradients(WaveFunction* waveFunction);
    Eigen::MatrixXd getAllInstantGradients();
    Eigen::MatrixXd updateParameters();
    Eigen::MatrixXd getEnergyGradient();

private:
    double m_A     = 20;
    double m_gamma = 0;
    double m_t     = m_A;
    int    m_step  = 0;
    int    m_numberOfBatches = 10;
    Eigen::MatrixXd m_v;
};
