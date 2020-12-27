#include "randomnormal.h"
#include "../system.h"

RandomNormalWeights::RandomNormalWeights(System *system, const double variance)
    : InitialWeights(system)
{
    m_system = system;
    m_variance = variance;
}

void RandomNormalWeights::setupInitialWeights()
{
    m_numberOfElements = m_system->getNumberOfElements();
    m_maxParameters = m_system->getMaxParameters();
    m_parameters = 0.1 * m_system->getRandomNumberGenerator()->randomNormalMatrix(m_numberOfElements, m_maxParameters, 0, m_variance);
    m_system->updateAllParameters(m_parameters);
}

Eigen::MatrixXd RandomNormalWeights::getParameters()
{
    return m_parameters;
}
