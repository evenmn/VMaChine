#pragma once
#include "activation.h"

class LeakyReLU : public Activation
{
public:
    LeakyReLU(System *system, double a = 0.1);
    Eigen::VectorXd evaluate(Eigen::VectorXd x);
    Eigen::VectorXd gradient(Eigen::VectorXd x);
    Eigen::VectorXd laplacian(Eigen::VectorXd x);

private:
    double m_a = 0;
};
