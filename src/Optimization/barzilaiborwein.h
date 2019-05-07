#pragma once
#include "optimization.h"
#include <Eigen/Dense>

class BarzilaiBorwein : public Optimization {
public:
    BarzilaiBorwein(System* system);
    unsigned int    getNumberOfBatches() { return m_numberOfBatches; }
    Eigen::MatrixXd updateParameters();
    Eigen::MatrixXd getEnergyGradient();

private:
    unsigned int    m_iter            = 0;
    unsigned int    m_numberOfBatches = 1;

    double          m_omega           = 0;
    double          m_omega_sqrd      = 0;

    Eigen::MatrixXd m_oldParameters;
    Eigen::MatrixXd m_parameters;
    Eigen::MatrixXd m_oldGradients;
    Eigen::MatrixXd m_gradients;
};
