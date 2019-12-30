#include "dense.h"
#include "../Activation/activation.h"

Dense::Dense(System *system, int numberOfUnits, Activation *activation)
    : Layer(system)
{
    m_system = system;
    m_numberOfUnits = numberOfUnits;
    m_activation = activation;
}

void Dense::updateWeights(Eigen::MatrixXd WNew)
{
    m_W = WNew;
}

void Dense::initialize(int numberOfUnitsInPreviousLayer)
{
    m_numberOfUnitsInPreviousLayer = numberOfUnitsInPreviousLayer;
    m_W = Eigen::MatrixXd::Zero(numberOfUnitsInPreviousLayer + 1, m_numberOfUnits); // Add bias weights

    m_a = Eigen::VectorXd::Ones(m_numberOfUnits + 1);
    m_da = Eigen::VectorXd::Ones(m_numberOfUnits + 1);
    m_dda = Eigen::VectorXd::Ones(m_numberOfUnits + 1);
}

Layer::Vector2l Dense::getWeightDim() {
    Vector2l size;
    size(0) = m_W.rows();
    size(1) = m_W.cols();
    return size;
}

Eigen::VectorXd Dense::evaluate(Eigen::VectorXd a0)
{
    m_a0 = a0;
    m_z = a0.transpose() * m_W;
    return m_z;
}

Eigen::VectorXd Dense::activate()
{
    m_a.tail(m_numberOfUnits) = m_activation->evaluate(m_z);
    return m_a;
}

Eigen::VectorXd Dense::activateDer() {
    m_da.tail(m_numberOfUnits) = m_activation->gradient(m_z);
    return m_da;
}

Eigen::VectorXd Dense::activateSecDer() {
    m_dda.tail(m_numberOfUnits) = m_activation->laplacian(m_z);
    return m_dda;
}

Eigen::VectorXd Dense::calculateDelta(Eigen::VectorXd delta1)
{
    Eigen::VectorXd da = m_da.tail(m_numberOfUnits);
    m_delta = m_W * da.cwiseProduct(delta1).tail(m_numberOfUnitsInPreviousLayer);
    return m_delta;
}

Eigen::MatrixXd Dense::calculateParameterGradient()
{
    Eigen::VectorXd da = m_da.tail(m_numberOfUnits);
    return m_delta.cwiseProduct(da) * m_a0.transpose();
}
