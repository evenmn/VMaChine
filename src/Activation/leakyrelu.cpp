#include "leakyrelu.h"
#include "../system.h"
#include <iostream>

LeakyReLU::LeakyReLU(System *system, double a)
    : Activation(system)
{
    m_system = system;
    m_a = a;
}

Eigen::VectorXd LeakyReLU::evaluate(Eigen::VectorXd x)
{
    return (x.array() > 0).select(x, m_a * x);
}

Eigen::VectorXd LeakyReLU::gradient(Eigen::VectorXd x)
{
    Eigen::VectorXd zeros = Eigen::VectorXd::Zero(x.size());
    Eigen::VectorXd ones = Eigen::VectorXd::Ones(x.size());
    return (x.array() > 0).select(m_a * ones, zeros);
}

Eigen::VectorXd LeakyReLU::laplacian(Eigen::VectorXd x)
{
    return Eigen::VectorXd::Zero(x.size());
}
