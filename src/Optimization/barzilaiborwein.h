#pragma once
#include "optimization.h"

class BarzilaiBorwein : public Optimization {
public:
    BarzilaiBorwein(System* system);
    int             getNumberOfBatches() { return m_numberOfBatches; }
    std::string     getLabel()           { return m_label; }
    Eigen::MatrixXd updateParameters();

private:
    int m_step, m_numberOfBatches, m_numberOfFreeDimensions;
    double m_omega, m_omega_sqrd;
    Eigen::MatrixXd m_oldParameters;
    Eigen::MatrixXd m_parameters;
    Eigen::MatrixXd m_oldGradients;
    Eigen::MatrixXd m_gradients;

    std::string     m_label = "BB";
};
