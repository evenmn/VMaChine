#include "input.h"
#include "../Activation/activation.h"

Input::Input(System *system, int numberOfUnits)
    : Layer(system)
{
    m_system = system;
    m_numberOfUnits = numberOfUnits;
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

void Input::initialize(int /*numberOfUnitsInPreviousLayer*/)
{
    m_a = Eigen::VectorXd::Ones(m_numberOfUnits + 1);
}

Eigen::VectorXd Input::evaluate(Eigen::VectorXd x)
{
    m_x = x;
    return m_x;
}

Eigen::VectorXd Input::activate()
{
    m_a.tail(m_numberOfUnits) = m_x;
    return m_a;
}

Eigen::VectorXd Input::activateDer()
{
    return Eigen::VectorXd::Zero(m_numberOfUnits + 1);
}

Eigen::VectorXd Input::activateSecDer()
{
    return Eigen::VectorXd::Zero(m_numberOfUnits + 1);
}

Eigen::VectorXd Input::calculateDelta(Eigen::VectorXd delta1)
{
    Eigen::VectorXd df = m_activation->gradient(m_x);
    m_delta = df.cwiseProduct(m_W * delta1);
    return m_delta;
}
