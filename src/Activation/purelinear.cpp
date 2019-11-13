#include "purelinear.h"
#include "../system.h"
#include <iostream>

PureLinear::PureLinear(System *system)
    : Activation(system)
{
    m_system = system;
}

Eigen::VectorXd PureLinear::evaluate(Eigen::VectorXd x)
{
    return x;
}

Eigen::VectorXd PureLinear::gradient(Eigen::VectorXd x)
{
    return Eigen::VectorXd::Ones(x.size());
}

Eigen::VectorXd PureLinear::laplacian(Eigen::VectorXd x)
{
    return Eigen::VectorXd::Zero(x.size());
}
