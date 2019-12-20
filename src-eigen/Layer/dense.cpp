#include "dense.h"
#include "../Activation/activation.h"

Dense::Dense(System *system, int numberOfUnits, Activation *activation)
    : Layer(system)
{
    m_system = system;
    m_numberOfUnits = numberOfUnits;
    m_activation = activation;
    //initialize();
}

void Dense::updateWeights(Eigen::MatrixXd WNew)
{
    m_W = WNew;
}

void Dense::initialize(int numberOfUnitsInPreviousLayer)
{
    m_numberOfUnitsInPreviousLayer = numberOfUnitsInPreviousLayer;
    m_W = Eigen::MatrixXd::Random(numberOfUnitsInPreviousLayer + 1, m_numberOfUnits); // Add bias weights
}

Layer::Vector2l Dense::getWeightDim() {
    Vector2l size;
    size(0) = m_W.rows();
    size(1) = m_W.cols();
    return size;
}

Eigen::VectorXd Dense::evaluate()
{
    m_z = m_a0.transpose() * m_W;
    return m_z;
}

Eigen::VectorXd Dense::activate(Eigen::VectorXd a0)
{
    m_a0 = Eigen::VectorXd::Ones(m_numberOfUnitsInPreviousLayer + 1); // Add bias unit
    m_a0.tail(m_numberOfUnitsInPreviousLayer) = a0;
    evaluate();
    m_a = m_activation->evaluate(m_z);
    return m_a;
}

Eigen::VectorXd Dense::calculateDelta(Eigen::VectorXd delta1)
{
    Eigen::VectorXd df = m_activation->gradient(m_z);
    m_delta = df.cwiseProduct(m_W * delta1);
    return m_delta;
}

Eigen::MatrixXd Dense::calculateGradient()
{
    Eigen::MatrixXd dC = m_a0 * m_delta.transpose();
    return dC;
}
