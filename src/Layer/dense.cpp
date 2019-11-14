#include "dense.h"
#include "../Activation/activation.h"

Dense::Dense(System *system, int h0, int h1, Activation *activation)
    : Layer(system)
{
    m_system = system;
    m_h0 = h0;
    m_h1 = h1;
    m_activation = activation;
    //initialize();
}

void Dense::updateWeights(Eigen::MatrixXd WNew)
{
    m_W = WNew;
}

void Dense::initialize()
{
    m_W = Eigen::MatrixXd::Random(m_h0 + 1, m_h1); // Add bias weights
}

Vector2l Dense::getWeightDim() {
    Vector2l size;
    size(0) = m_W.rows();
    size(1) = m_W.cols();
    return size;
}

Eigen::VectorXd Dense::evaluate(Eigen::VectorXd a0)
{
    Eigen::VectorXd a = Eigen::VectorXd::Ones(m_h0 + 1); // Add bias unit
    a.tail(m_h0) = a0;
    Eigen::VectorXd z = a * m_W;
    return z;
}

Eigen::VectorXd Dense::activate(Eigen::VectorXd a0)
{
    Eigen::VectorXd z = evaluate(a0);
    Eigen::VectorXd a = m_activation->evaluate(z);
    return a;
}

Eigen::VectorXd Dense::activateDer(Eigen::VectorXd z)
{
    Eigen::VectorXd df = m_activation->gradient(z);
    return df;
}

Eigen::VectorXd Dense::activateSecDer(Eigen::VectorXd z)
{
    Eigen::VectorXd ddf = m_activation->laplacian(z);
    return ddf;
}

Eigen::VectorXd Dense::calculateDelta(Eigen::VectorXd z)
{
    Eigen::VectorXd df = activateDer(z);
    Eigen::VectorXd delta = df;
    return delta;
}

Eigen::MatrixXd Dense::calculateGradient(Eigen::VectorXd z, Eigen::VectorXd a0)
{
    Eigen::VectorXd delta = calculateDelta(z);
    Eigen::MatrixXd dC = a0.transpose() * delta;
    return dC;
}
