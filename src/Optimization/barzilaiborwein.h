#pragma once
#include "optimization.h"

class BarzilaiBorwein : public Optimization {
public:
    BarzilaiBorwein(System* system);
    int             getNumberOfBatches() { return m_numberOfBatches; }
    std::string     getLabel()           { return m_label; }
    Eigen::MatrixXd updateParameters();

private:
    int    m_step = 0;
    int    m_numberOfBatches = 1;
    double m_omega = 0;
    double m_omega_sqrd = 0;
    Eigen::MatrixXd m_oldParameters;
    Eigen::MatrixXd m_parameters;
    Eigen::MatrixXd m_oldGradients;
    Eigen::MatrixXd m_gradients;

    std::string     m_label = "BB";
};
