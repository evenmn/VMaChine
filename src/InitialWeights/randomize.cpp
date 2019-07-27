#include "randomize.h"
#include <iostream>
#include <cassert>
#include "../system.h"

Randomize::Randomize(System* system, const double factor = 1) :
    InitialWeights(system) {
    m_system                = system;
    m_numberOfDimensions    = m_system->getNumberOfDimensions();
    m_numberOfParticles     = m_system->getNumberOfParticles();
    m_numberOfElements      = m_system->getNumberOfElements();
    m_maxParameters         = m_system->getMaxParameters();
    m_factor                = factor;
    setupInitialWeights();
}

void Randomize::setupInitialWeights() {
    m_parameters = m_factor * Eigen::MatrixXd::Random(m_numberOfElements, m_maxParameters);
    m_parameters(3,0) = 1;
    m_system->updateAllParameters(m_parameters);
}

Eigen::MatrixXd Randomize::getParameters() {
    return m_parameters;
}
