#include "constant.h"
#include <iostream>
#include <cassert>
#include "../system.h"

Constant::Constant(System* system, const double factor)  :  InitialWeights(system) {
    m_system                          = system;
    m_numberOfDimensions              = m_system->getNumberOfDimensions();
    m_numberOfParticles               = m_system->getNumberOfParticles();
    m_numberOfElements                = m_system->getNumberOfElements();
    m_maxNumberOfParametersPerElement = m_system->getMaxParameters();
    m_factor                          = factor;
    setupInitialWeights();
}

void Constant::setupInitialWeights() {
    m_parameters = m_factor * Eigen::MatrixXd::Ones(Eigen::Index(m_numberOfElements), m_maxNumberOfParametersPerElement);
    m_system->updateAllParameters(m_parameters);
}
