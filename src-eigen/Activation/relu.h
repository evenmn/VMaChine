#pragma once
#include "activation.h"

class ReLU : public Activation
{
public:
    ReLU(System *system, double a = 0);
    Eigen::VectorXd evaluate(Eigen::VectorXd x);
    Eigen::VectorXd gradient(Eigen::VectorXd x);
    Eigen::VectorXd laplacian(Eigen::VectorXd x);

private:
    double m_a = 0;
};
