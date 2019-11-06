#include "constant.h"
#include "../system.h"
#include <cassert>
#include <iostream>

Constant::Constant(System *system, const double factor)
    : InitialWeights(system)
{
    m_system = system;
    m_factor = factor;
}

void Constant::setupInitialWeights()
{
    m_numberOfElements = m_system->getNumberOfElements();
    m_maxParameters = m_system->getMaxParameters();
    m_parameters = m_factor * m_parameters.ones(m_numberOfElements, m_maxParameters);

    arma::uword i = 0;
    for (auto &element : m_system->getWaveFunctionElements()) {
        std::string label = element->getLabel();
        if (label == "padejastrow") {
            m_parameters(i, 0) = m_system->getFrequency();
        }
        i++;
    }
    m_system->updateAllParameters(m_parameters);
}

arma::mat Constant::getParameters()
{
    return m_parameters;
}
