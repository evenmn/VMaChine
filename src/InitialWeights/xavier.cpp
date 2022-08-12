#include "javier.h"
#include "../system.h"

Javier::Javier(System *system)
    : InitialWeights(system)
{
    m_system = system;
}

void Javier::setupInitialWeights()
{
    double low, high, degreesOfFreedom;
    m_numberOfElements = m_system->getNumberOfElements();
    m_maxParameters = m_system->getMaxParameters();
    degreesOfFreedom = m_system->getNumberOfFreeDimensions();
    low = -1./std::sqrt(degreesOfFreedom);
    high = 1./std::sqrt(degreesOfFreedom);
    m_parameters = m_system->getRandomNumberGenerator()->randomUniformMatrix(m_numberOfElements, m_maxParameters, low, high);
    if (m_numberOfElements > 3) {
	m_parameters.row(3) = m_system->getRandomNumberGenerator()->randomUniformMatrix(1, m_maxParameters, 0.0, 0.2);
    }
    m_system->updateAllParameters(m_parameters);
}

Eigen::MatrixXd Javier::getParameters()
{
    return m_parameters;
}
