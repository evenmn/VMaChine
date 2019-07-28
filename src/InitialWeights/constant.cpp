#include "constant.h"
#include <iostream>
#include <cassert>
#include "../system.h"

Constant::Constant(System* system, const double factor = 1) :
    InitialWeights(system) {
    m_system                = system;
    m_numberOfDimensions    = m_system->getNumberOfDimensions();
    m_numberOfParticles     = m_system->getNumberOfParticles();
    m_numberOfElements      = m_system->getNumberOfElements();
    m_maxParameters         = m_system->getMaxParameters();
    m_factor                = factor;
    setupInitialWeights();
}

void Constant::setupInitialWeights() {
    m_parameters = m_factor * Eigen::MatrixXd::Ones(m_numberOfElements, m_maxParameters);
    int i = 0;
    for(auto& element : m_system->getWaveFunctionElements()) {
        std::string label = element->getLabel();
        if(label == "padejastrow") {
            m_parameters(i,0) = 1;
        }
        i++;
    }
    m_system->updateAllParameters(m_parameters);
}

Eigen::MatrixXd Constant::getParameters() {
    return m_parameters;
}
