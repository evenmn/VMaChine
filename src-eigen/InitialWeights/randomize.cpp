#include "randomize.h"
#include "../system.h"
#include <cassert>
#include <iostream>

Randomize::Randomize(System *system, const double factor)
    : InitialWeights(system)
{
    m_system = system;
    m_factor = factor;
}

void Randomize::setupInitialWeights()
{
    m_numberOfElements = m_system->getNumberOfElements();
    m_maxParameters = m_system->getMaxParameters();
    m_parameters = m_factor * Eigen::MatrixXd::Random(m_numberOfElements, m_maxParameters);
    m_system->updateAllParameters(m_parameters);
}

Eigen::MatrixXd Randomize::getParameters()
{
    return m_parameters;
}
