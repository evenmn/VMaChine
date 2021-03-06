#pragma once
#include "optimization.h"

class ADAM : public Optimization
{
public:
    ADAM(System *system);
    int getNumberOfBatches() { return m_numberOfBatches; }
    std::string getLabel() { return m_label; }

    void initialize();
    Eigen::MatrixXd updateParameters();

private:
    int m_numberOfBatches = 10;
    double m_beta1 = 0.89;
    double m_beta2 = 0.999;
    double m_epsilon = 1e-8;
    int m_step = 0;
    Eigen::MatrixXd m_g;
    Eigen::MatrixXd m_m;
    Eigen::MatrixXd m_v;
    Eigen::MatrixXd m_mHat;
    Eigen::MatrixXd m_vHat;
    Eigen::MatrixXd m_theta;
    Eigen::MatrixXd m_denom;

    std::string m_label = "ADAM";
};
