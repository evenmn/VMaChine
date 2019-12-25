#include "input.h"
#include "../Activation/activation.h"

Input::Input(System *system, int numberOfUnits)
    : Layer(system)
{
    m_system = system;
    m_numberOfUnits = numberOfUnits;
    //initialize();
}

void Input::updateWeights(Eigen::MatrixXd WNew)
{
    m_W = WNew;
}

Layer::Vector2l Input::getWeightDim()
{
    Vector2l size;
    size(0) = m_W.rows();
    size(1) = m_W.cols();
    return size;
}

void Input::initialize(int numberOfUnitsInPreviousLayer, double factor)
{
    int g = numberOfUnitsInPreviousLayer;
}

Eigen::VectorXd Input::evaluate(Eigen::VectorXd a0)
{
    m_x = a0;
    return m_x;
}

Eigen::VectorXd Input::activate()
{
    return m_x;
}

Eigen::VectorXd Input::activateDer()
{
    return m_x;
}

Eigen::VectorXd Input::activateSecDer()
{
    return m_x;
}

Eigen::VectorXd Input::calculateDelta(Eigen::VectorXd delta1)
{
    Eigen::VectorXd df = m_activation->gradient(m_x);
    m_delta = df.cwiseProduct(m_W * delta1);
    return m_delta;
}

Eigen::MatrixXd Input::calculateGradient()
{
    Eigen::MatrixXd dC = m_x * m_delta.transpose();
    return dC;
}
