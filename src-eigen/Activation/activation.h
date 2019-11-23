#pragma once
#include "../Eigen/Dense"
#include <iostream>

class Activation
{
public:
    Activation(class System *system);
    virtual Eigen::VectorXd evaluate(Eigen::VectorXd x) = 0;
    virtual Eigen::VectorXd gradient(Eigen::VectorXd x) = 0;
    virtual Eigen::VectorXd laplacian(Eigen::VectorXd x) = 0;

    virtual ~Activation() = 0;

protected:
    class System *m_system = nullptr;
};
