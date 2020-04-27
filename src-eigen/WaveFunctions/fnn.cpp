#include "fnn.h"
#include "../system.h"

FNN::FNN(System *system)
    : WaveFunction(system)
{}

void FNN::setConstants(const int elementNumber)
{
    m_elementNumber = elementNumber;
    m_degreesOfFreedom = m_system->getNumberOfFreeDimensions();
    m_layers2 = m_system->getLayers();

    // Set up weight matrices
    int h0 = m_layers2[0]->getNumberOfUnits();
    m_layers2[0]->initialize(0);
    //int h1 = m_layers2[1]->getNumberOfUnits();
    //m_layers2[1]->initialize(h0);
    //int h2 = m_layers2[2]->getNumberOfUnits();
    //m_layers2[2]->initialize(h1);
    //m_numberOfParameters += (h0 + 1) * h1;
    //m_numberOfParameters += (h1 + 1) * h2;
    for(unsigned long l = 1; l < m_layers2.size(); l++) {
        m_layers2[l]->initialize(h0);
        int h1 = m_layers2[l]->getNumberOfUnits();
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
    //std::cout << "evaluate1" << std::endl;
    Eigen::VectorXd a = position;
    for(unsigned long l = 0; l < m_layers2.size(); l++) {
        m_layers2[l]->evaluate(a);
        a = m_layers2[l]->activate();
        m_layers2[l]->activateDer();
        m_layers2[l]->activateSecDer();
    }
    //std::cout << "evaluate2" << std::endl;
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
    //std::cout << "updateArrays1" << std::endl;
    m_positions = positions;
    m_out = evaluate(positions);
    updateProbabilityRatio(changedCoord);
    //std::cout << "updateArrays2" << std::endl;
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
    //std::cout << "updateParameters1" << std::endl;
    int cumulativeStart = 0;
    for(unsigned long i = 0; i < m_layers2.size()-1; i++) {
        Layer::Vector2l size = m_layers2[i+1]->getWeightDim();
        Eigen::VectorXd WFlatten = parameters.row(m_elementNumber)
                                       .segment(cumulativeStart, size.prod());
        Eigen::MatrixXd W = WaveFunction::reshape(WFlatten, size(0), size(1));
        m_layers2[i+1]->updateWeights(W);
        cumulativeStart += size.prod();
    }
    //std::cout << "updateParameters2" << std::endl;
}

double FNN::evaluateRatio()
{
    return m_probabilityRatio * m_probabilityRatio;
}

double FNN::computeGradient(const int k)
{
    //std::cout << "computeGradient1" << std::endl;
    //Eigen::VectorXd delta = Eigen::VectorXd::Ones(1);
    //for(unsigned long l = m_layers2.size()-1; l > 0; l--) {
    //    delta = m_layers2[l]->calculateDelta(delta);
    //}
    //return delta(k) / m_out;

    Eigen::VectorXd da1 = m_layers2[1]->getDA();
    Eigen::VectorXd da2 = m_layers2[2]->getDA();
    Eigen::MatrixXd w1 = m_layers2[1]->getWeights();
    Eigen::MatrixXd w2 = m_layers2[2]->getWeights();

    double gradient = 0;
    for(int i=0; i < da1.size() - 1; i++) {
        for(int j=0; j < da2.size() - 1; j++){
            gradient += da2(j) * da1(i) * w2(i, j) * w1(k, i);
        }
    }
    //return gradient / m_out;
    //std::cout << "computeGradient2" << std::endl;
    return 0;
}

double FNN::computeLaplacian()
{
    //std::cout << "computeLaplacian1" << std::endl;
    Eigen::VectorXd da1 = m_layers2[1]->getDA();
    Eigen::VectorXd dda1 = m_layers2[1]->getDDA();
    Eigen::VectorXd da2 = m_layers2[2]->getDA();
    Eigen::VectorXd dda2 = m_layers2[2]->getDDA();
    Eigen::MatrixXd w1 = m_layers2[1]->getWeights();
    Eigen::MatrixXd w2 = m_layers2[2]->getWeights();

    double laplacian = 0;
    double gradientSqrd = 0;
    for(int k=0; k < m_degreesOfFreedom; k++) {
        for(int j=0; j < dda2.size() - 1; j++) {
            for(int h=0; h < da1.size() - 1; h++) {
                laplacian += dda2(j) * da1(h) * w1(k, h) * w2(h, j) * da1(h) * w1(k, h) * w2(h, j);
                laplacian += da2(j) * dda1(h) * w1(k, h) * w2(h, j) * w1(k, h);
            }
        }
        double compGradient = computeGradient(k);
        gradientSqrd += compGradient * compGradient;
    }
    //std::cout << laplacian << std::endl;
    //std::cout << gradientSqrd << std::endl;
    //std::cout << std::endl;
    //std::cout << "computeLaplacian2" << std::endl;
    return laplacian / m_out - gradientSqrd / (m_out * m_out);
}

Eigen::VectorXd FNN::computeParameterGradient()
{
    //std::cout << "computeParameterGradient1" << std::endl;
    m_gradients = Eigen::VectorXd::Zero(m_system->getMaxParameters());



//    Eigen::VectorXd delta = Eigen::VectorXd::Ones(1);
//    int cumulativeStart = 0;
//    for (unsigned long l = m_layers2.size()-1; l > 0; l--) {
//        Layer::Vector2l size = m_layers2[l]->getWeightDim();
//        Eigen::VectorXd da = m_layers2[l]->getDA().tail(size(1));
//        Eigen::VectorXd a0 = m_layers2[l-1]->getA();
//        Eigen::MatrixXd dW = delta.cwiseProduct(da) * a0.transpose() / m_out;
//        std::cout << dW << std::endl;
//        m_gradients.segment(cumulativeStart, size.prod()) = WaveFunction::flatten(dW);
//        std::cout << WaveFunction::flatten(dW) << std::endl;
//        std::cout << "\n\n\n\n\n" << std::endl;
//        delta = m_layers2[l]->calculateDelta(delta);
//        cumulativeStart += size.prod();
//    }
    m_gradients = Eigen::VectorXd::Zero(m_system->getMaxParameters());
    //std::cout << "computeParameterGradient2" << std::endl;
    return m_gradients;
}
