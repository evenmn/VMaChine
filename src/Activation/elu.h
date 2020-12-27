#pragma once
#include "activation.h"

class ELU : public Activation
{
public:
    ELU(System *system);
    Eigen::VectorXd evaluate(Eigen::VectorXd x);
    Eigen::VectorXd gradient(Eigen::VectorXd x);
    Eigen::VectorXd laplacian(Eigen::VectorXd x);
};
