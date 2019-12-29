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
    m_out = evaluate(positions);
}

double FNN::evaluate(Eigen::VectorXd position) {
    Eigen::VectorXd a = position;
    for(unsigned long i = 1; i < m_layers2.size(); i++) {
        m_layers2[i]->evaluate(a);
        a = m_layers2[i]->activate();
        //std::cout << a << std::endl;
        //std::cout << std::endl;
        m_layers2[i]->activateDer();
        m_layers2[i]->activateSecDer();
    }
    return a(0);
}

void FNN::updateProbabilityRatio(int /*changedCoord*/)
{   
    m_probabilityRatio = m_out / m_outOld;
}

void FNN::updateArrays(const Eigen::VectorXd positions,
                       const Eigen::VectorXd /*radialVector*/,
                       const Eigen::MatrixXd /*distanceMatrix*/,
                       const int changedCoord)
{
    m_positions = positions;
    updateProbabilityRatio(changedCoord);
    evaluate(positions);
}

void FNN::setArrays()
{
    m_positionsOld = m_positions;
    m_probabilityRatioOld = m_probabilityRatio;
    m_outOld = m_out;
}

void FNN::resetArrays()
{
    m_positions = m_positionsOld;
    m_probabilityRatio = m_probabilityRatioOld;
    m_out = m_outOld;
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
    Eigen::VectorXd delta = Eigen::VectorXd::Ones(1);
    for(unsigned long l = m_layers2.size()-1; l > 0; l--) {
        delta = m_layers2[l]->calculateDelta(delta);
    }
    return delta(k) / m_out;
}

double FNN::computeLaplacian()
{
    Eigen::VectorXd da1 = m_layers2[1]->getDA();
    Eigen::VectorXd dda1 = m_layers2[1]->getDDA();
    Eigen::VectorXd da2 = m_layers2[2]->getDA();
    Eigen::VectorXd dda2 = m_layers2[2]->getDDA();
    Eigen::MatrixXd w1 = m_layers2[1]->getWeights();
    Eigen::MatrixXd w2 = m_layers2[2]->getWeights();

    double result = 0;
    for(int k=0; k < m_degreesOfFreedom; k++) {
        for(int j=0; j < dda2.size(); j++) {
            for(int h=0; h < da1.size(); h++) {
                result += dda2(j) * da1(h) * w1(k, h) * w2(h, j) * da1(h) * w1(k, h) * w2(h, j);
                result += da2(j) * dda1(h) * w1(k, h) * w2(h, j) * w1(k, h) * w2(h, j);
            }
        }
    }
    return result / m_out;
}

Eigen::VectorXd FNN::computeParameterGradient()
{
    m_gradients = Eigen::VectorXd::Zero(m_system->getMaxParameters());
    Eigen::VectorXd da2 = m_layers2[2]->getDA();
    Eigen::VectorXd a1 = m_layers2[1]->getA();

    m_gradients.head(5) = a1 * da2.transpose() / m_out;


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
