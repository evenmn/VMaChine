#include "randomuniform.h"
#include "../system.h"

RandomUniformWeights::RandomUniformWeights(System *system, const double factor)
    : InitialWeights(system)
{
    m_system = system;
    m_factor = factor;
}

void RandomUniformWeights::setupInitialWeights()
{
    m_numberOfElements = m_system->getNumberOfElements();
    m_maxParameters = m_system->getMaxParameters();
    m_parameters = m_factor * m_system->getRandomNumberGenerator()->randomUniformMatrix(m_numberOfElements, m_maxParameters);
    m_system->updateAllParameters(m_parameters);
}

Eigen::MatrixXd RandomUniformWeights::getParameters()
{
    return m_parameters;
}
