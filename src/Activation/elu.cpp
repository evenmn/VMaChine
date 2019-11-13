#include "elu.h"
#include "../system.h"
#include <iostream>

ELU::ELU(System *system)
    : Activation(system)
{
    m_system = system;
}

Eigen::VectorXd ELU::evaluate(Eigen::VectorXd x)
{
    return (x.array() > 0).select(x, x.array().exp() - 1);
}

Eigen::VectorXd ELU::gradient(Eigen::VectorXd x)
{
    Eigen::VectorXd ones = Eigen::VectorXd::Ones(x.size());
    return (x.array() > 0).select(ones, x.array().exp());
}

Eigen::VectorXd ELU::laplacian(Eigen::VectorXd x)
{
    Eigen::VectorXd zeros = Eigen::VectorXd::Zero(x.size());
    return (x.array() > 0).select(zeros, x.array().exp());
}
