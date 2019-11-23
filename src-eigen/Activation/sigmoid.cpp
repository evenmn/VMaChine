#include "sigmoid.h"
#include "../system.h"
#include <iostream>

Sigmoid::Sigmoid(System *system)
    : Activation(system)
{
    /* Constructor for the Sigmoid class. */
    m_system = system;
}

Eigen::VectorXd Sigmoid::evaluate(Eigen::VectorXd x)
{
    /* Evaluate the sigmoid activation function at
     * point x. */
    return 1 / (1 + x.array().exp());
}

Eigen::VectorXd Sigmoid::gradient(Eigen::VectorXd x)
{
    /* Evaluate the derivative (gradient) of the sigmoid
     * activation function at point x. */
    Eigen::VectorXd ones = Eigen::VectorXd::Ones(x.size());
    return evaluate(x).cwiseProduct(ones - evaluate(x));
}

Eigen::VectorXd Sigmoid::laplacian(Eigen::VectorXd x)
{
    /* Evaluate the second derivative (Laplacian) of the
     * sigmoid activation function at point x. */
    Eigen::VectorXd ones = Eigen::VectorXd::Ones(x.size());
    return gradient(x).cwiseProduct(ones - 2 * evaluate(x));
}
