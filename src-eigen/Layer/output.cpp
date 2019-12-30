#include "output.h"
#include "../Activation/activation.h"

Output::Output(System *system, Activation *activation)
    : Layer(system)
{
    m_system = system;
    m_numberOfUnits = 1;
    m_activation = activation;
}

void Output::updateWeights(Eigen::MatrixXd WNew)
{
    m_W = WNew;
}

void Output::initialize(int numberOfUnitsInPreviousLayer)
{
    m_numberOfUnitsInPreviousLayer = numberOfUnitsInPreviousLayer;
    m_W = Eigen::MatrixXd::Zero(numberOfUnitsInPreviousLayer + 1, m_numberOfUnits); // Add bias weights
}

Layer::Vector2l Output::getWeightDim() {
    Vector2l size;
    size(0) = m_W.rows();
    size(1) = m_W.cols();
    return size;
}

Eigen::VectorXd Output::evaluate(Eigen::VectorXd a0)
{
    m_z = a0.transpose() * m_W;
    return m_z;
}

Eigen::VectorXd Output::activate()
{
    m_a = m_activation->evaluate(m_z);
    return m_a;
}

Eigen::VectorXd Output::activateDer() {
    m_da = m_activation->gradient(m_z);
    return m_da;
}

Eigen::VectorXd Output::activateSecDer() {
    m_dda = m_activation->laplacian(m_z);
    return m_dda;
}

Eigen::VectorXd Output::calculateDelta(Eigen::VectorXd delta1)
{
    Eigen::VectorXd da = m_da.tail(m_numberOfUnits);
    m_delta = m_W * da.cwiseProduct(delta1);
    return m_delta.tail(m_numberOfUnitsInPreviousLayer);
}
