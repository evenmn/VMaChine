#pragma once
#include "optimization.h"
#include <Eigen/Dense>

class ASGD : public Optimization {
public:
    ASGD(System* system, const double gamma);
    unsigned int    getNumberOfBatches() { return m_numberOfBatches; }
    std::string     getLabel()           { return m_label; }
    Eigen::MatrixXd updateParameters();

private:
    unsigned int    m_iter            = 0;
    unsigned int    m_numberOfBatches = 10;

    double          m_A               = 20;
    double          m_gamma           = 0;
    double          m_t               = m_A;

    Eigen::MatrixXd m_v;

    std::string     m_label = "ASGD";
};
