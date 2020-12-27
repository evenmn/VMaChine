#pragma once
#include "activation.h"

class Sigmoid : public Activation
{
public:
    Sigmoid(System *system);
    Eigen::VectorXd evaluate(Eigen::VectorXd x);
    Eigen::VectorXd gradient(Eigen::VectorXd x);
    Eigen::VectorXd laplacian(Eigen::VectorXd x);
};
