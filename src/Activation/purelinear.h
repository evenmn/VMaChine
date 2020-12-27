#pragma once
#include "activation.h"

class PureLinear : public Activation
{
public:
    PureLinear(System *system);
    Eigen::VectorXd evaluate(Eigen::VectorXd x);
    Eigen::VectorXd gradient(Eigen::VectorXd x);
    Eigen::VectorXd laplacian(Eigen::VectorXd x);
};
