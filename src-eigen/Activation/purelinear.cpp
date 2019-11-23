#include "purelinear.h"
#include "../system.h"
#include <iostream>

PureLinear::PureLinear(System *system)
    : Activation(system)
{
    /* Constructor for the PureLinear class. */
    m_system = system;
}

Eigen::VectorXd PureLinear::evaluate(Eigen::VectorXd x)
{
    /* Evaluate the pure linear activation function
     * at point x. */
    return x;
}

Eigen::VectorXd PureLinear::gradient(Eigen::VectorXd x)
{
    /* Evaluate the derivative (gradient) of the
     * pure linear activation function at point x. */
    return Eigen::VectorXd::Ones(x.size());
}

Eigen::VectorXd PureLinear::laplacian(Eigen::VectorXd x)
{
    /* Evaluate the second derivative (Laplacian) of the
     * pure linear activation function at point x. */
    return Eigen::VectorXd::Zero(x.size());
}
