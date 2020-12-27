#include "relu.h"
#include "../system.h"
#include <iostream>

ReLU::ReLU(System *system, double a)
    : Activation(system)
{
    /* Constructor for the ReLU class. */
    m_system = system;
    m_a = a;
}

Eigen::VectorXd ReLU::evaluate(Eigen::VectorXd x)
{
    /* Evaluate the ReLU activation function at
     * point x. */
    return (x.array() > 0).select(x, m_a * x);
}

Eigen::VectorXd ReLU::gradient(Eigen::VectorXd x)
{
    /* Evaluate the derivative (gradient) of the
     * ReLU activation function at point x. */
    Eigen::VectorXd ones = Eigen::VectorXd::Ones(x.size());
    return (x.array() > 0).select(ones, m_a * ones);
}

Eigen::VectorXd ReLU::laplacian(Eigen::VectorXd x)
{
    /* Evaluate the second derivative (Laplacian) of the
     * ReLU activation function at point x. */
    return Eigen::VectorXd::Zero(x.size());
}
