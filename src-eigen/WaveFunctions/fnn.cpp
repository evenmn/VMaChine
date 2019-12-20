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
    m_layers2 = m_system->getLayers();

    // Set up weight matrices
    int h0 = m_layers2[0]->getNumberOfUnits();
    for(unsigned long i = 0; i < m_layers2.size()-1; i++) {
        m_layers2[i+1]->initialize(h0);
        int h1 = m_layers2[i+1]->getNumberOfUnits();
        m_numberOfParameters += (h0 + 1) * h1;
        h0 = h1;
    }
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
    for(unsigned long i = 0; i < m_layers2.size()-1; i++) {
        a = m_layers2[i+1]->activate(a);
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
    for(unsigned long i = 0; i < m_layers2.size()-1; i++) {
        Layer::Vector2l size = m_layers2[i+1]->getWeightDim();
        Eigen::VectorXd WFlatten = parameters.row(m_elementNumber)
                                       .segment(cumulativeStart, size.prod());
        Eigen::MatrixXd W = WaveFunction::reshape(WFlatten, size(0), size(1));
        m_layers2[i+1]->updateWeights(W);
        cumulativeStart += size.prod();
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
    /*
    int cumulativeStart = 0;
    for (unsigned long i = 0; i < m_layers2.size()-1; i++) {
        Layer::Vector2l size = m_layers2[i+1]->getWeightDim();
        std::cout << size.transpose() << std::endl;
        m_layers2[i+1]->calculateDelta()
        Eigen::MatrixXd W = m_layers2[i+1]->calculateGradient();
        std::cout << W << std::endl;
        m_gradients.segment(cumulativeStart, size.prod()) = WaveFunction::flatten(W);
        cumulativeStart += size.prod();
    }
    */
    return m_gradients;
}
