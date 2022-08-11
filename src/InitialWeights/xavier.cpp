#include "xavier.h"
#include "../system.h"

Xavier::Xavier(System *system)
    : InitialWeights(system)
{
    m_system = system;
}

void Xavier::setupInitialWeights()
{
    double low, high, degreesOfFreedom;
    m_numberOfElements = m_system->getNumberOfElements();
    m_maxParameters = m_system->getMaxParameters();
    degreesOfFreedom = m_system->getNumberOfFreeDimensions();
    low = -1./std::sqrt(degreesOfFreedom);
    high = 1./std::sqrt(degreesOfFreedom);
    m_parameters = m_system->getRandomNumberGenerator()->randomUniformMatrix(m_numberOfElements, m_maxParameters, low, high);
    m_system->updateAllParameters(m_parameters);
}

Eigen::MatrixXd Xavier::getParameters()
{
    return m_parameters;
}
