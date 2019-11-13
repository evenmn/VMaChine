#include "sigmoid.h"
#include "../system.h"
#include <iostream>

Sigmoid::Sigmoid(System *system)
    : Activation(system)
{
    m_system = system;
}

Eigen::VectorXd Sigmoid::evaluate(Eigen::VectorXd x)
{
    return 1 / (1 + x.array().exp());
}

Eigen::VectorXd Sigmoid::gradient(Eigen::VectorXd x)
{
    Eigen::VectorXd ones = Eigen::VectorXd::Ones(x.size());
    return evaluate(x).cwiseProduct(ones - evaluate(x));
}

Eigen::VectorXd Sigmoid::laplacian(Eigen::VectorXd x)
{
    Eigen::VectorXd ones = Eigen::VectorXd::Ones(x.size());
    return gradient(x).cwiseProduct(ones - 2 * evaluate(x));
}
