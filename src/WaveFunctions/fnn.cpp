#include "fnn.h"
#include "../system.h"
#include <cassert>
#include <iostream>

FNN::FNN(System *system)
    : WaveFunction(system)
{}

void FNN::setConstants(const int elementNumber)
{
    m_elementNumber = elementNumber;
    m_degreesOfFreedom = m_system->getNumberOfFreeDimensions();
    m_omega = m_system->getFrequency();
    m_layers = m_system->getLayers();
    m_activationFunctions = m_system->getActivationFunctions();
    m_units = m_system->getHiddenUnits();
    m_units.insert(m_units.begin(), m_degreesOfFreedom);    // Add input layer
    m_units.push_back(1);                                   // Add output layer
    m_activationFunctions.push_back(new Sigmoid(m_system)); // Add activation on output
    m_numberOfParameters = m_numberOfHiddenUnits * (1 + m_degreesOfFreedom);
}

void FNN::initializeArrays(const Eigen::VectorXd positions,
                           const Eigen::VectorXd /*radialVector*/,
                           const Eigen::MatrixXd /*distanceMatrix*/)
{
    m_positions = positions;
    m_probabilityRatio = 1;
}

double FNN::evaluate(Eigen::VectorXd position) {
    Eigen::VectorXd a = position;
    for(auto &i : m_layers) {
        a = i->activate(a);
    }
    return a(0);
}

void FNN::updateProbabilityRatio(int /*changedCoord*/)
{   
    m_probabilityRatio = evaluate(m_positions) / evaluate(m_positionsOld);
}

void FNN::updateArrays(const Eigen::VectorXd positions,
                       const Eigen::VectorXd /*radialVector*/,
                       const Eigen::MatrixXd /*distanceMatrix*/,
                       const int changedCoord)
{
    m_positions = positions;
    updateProbabilityRatio(changedCoord);
}

void FNN::setArrays()
{
    m_positionsOld = m_positions;
    m_probabilityRatioOld = m_probabilityRatio;
}

void FNN::resetArrays()
{
    m_positions = m_positionsOld;
    m_probabilityRatio = m_probabilityRatioOld;
}

void FNN::updateParameters(const Eigen::MatrixXd parameters)
{
    int cumulativeStart = 0;
    for(auto &i : m_layers) {
        Layer::Vector2l size = i->getWeightDim();
        cumulativeStart += size.prod();
        Eigen::MatrixXd WFlatten = parameters.row(m_elementNumber)
                                       .segment(cumulativeStart, size.prod());
        i->updateWeights(WaveFunction::reshape(WFlatten, size(0), size(1)));
    }
}

double FNN::evaluateRatio()
{
    return m_probabilityRatio * m_probabilityRatio;
}

double FNN::computeGradient(const int k)
{
    return -m_omegalpha * m_positions(k);
}

double FNN::computeLaplacian()
{
    return -m_omegalpha * m_degreesOfFreedom;
}

Eigen::VectorXd FNN::computeParameterGradient()
{
    m_gradients = Eigen::VectorXd::Zero(m_system->getMaxParameters());

    int cumulativeStart = 0;
    for (auto &i : m_layers) {
        Layer::Vector2l size = i->getWeightDim();
        cumulativeStart += size.prod();
        Eigen::MatrixXd W = i->calculateGradient();
        m_gradients.segment(cumulativeStart, size.prod()) = WaveFunction::flatten(W);
    }
    return m_gradients;
}
