#include "relu.h"
#include "../system.h"
#include <iostream>

ReLU::ReLU(System *system, double a)
    : Activation(system)
{
    m_system = system;
    m_a = a;
}

Eigen::VectorXd ReLU::evaluate(Eigen::VectorXd x)
{
    return (x.array() > 0).select(x, m_a * x);
}

Eigen::VectorXd ReLU::gradient(Eigen::VectorXd x)
{
    Eigen::VectorXd zeros = Eigen::VectorXd::Zero(x.size());
    Eigen::VectorXd ones = Eigen::VectorXd::Ones(x.size());
    return (x.array() > 0).select(m_a * ones, zeros);
}

Eigen::VectorXd ReLU::laplacian(Eigen::VectorXd x)
{
    return Eigen::VectorXd::Zero(x.size());
}
