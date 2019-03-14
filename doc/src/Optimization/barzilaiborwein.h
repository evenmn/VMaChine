#pragma once
#include "optimization.h"
#include <Eigen/Dense>

class BarzilaiBorwein : public Optimization {
public:
    BarzilaiBorwein(System* system);
    int             getNumberOfBatches() { return m_numberOfBatches; }
    Eigen::VectorXd getInstantGradients(WaveFunction* waveFunction);
    Eigen::MatrixXd getAllInstantGradients();
    Eigen::MatrixXd updateParameters();
    Eigen::MatrixXd getEnergyGradient();

private:
    int    m_step = 0;
    int    m_numberOfBatches = 1;
    double m_omega = 0;
    double m_omega_sqrd = 0;
    Eigen::MatrixXd m_oldParameters;
    Eigen::MatrixXd m_parameters;
    Eigen::MatrixXd m_oldGradients;
    Eigen::MatrixXd m_gradients;
};
