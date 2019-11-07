#include "sgd.h"
#include "../WaveFunctions/wavefunction.h"
#include "../sampler.h"
#include "../system.h"
#include <cassert>
#include <iostream>

SGD::SGD(System *system, const double gamma, const double monotonicExp)
    : Optimization(system)
{
    m_gamma = gamma;
    m_monotonicExp = monotonicExp;
}

void SGD::initialize()
{
    m_numberOfElements = m_system->getNumberOfElements();
    m_maxParameters = m_system->getMaxParameters();
    m_eta = m_system->getLearningRate();
    m_v.ones(m_numberOfElements, m_maxParameters);
}

arma::mat SGD::updateParameters()
{
    m_step += 1;
    double monotonic = 1 / pow(m_step, m_monotonicExp);
    m_v = m_gamma * m_v + m_eta * Optimization::getEnergyGradient() * monotonic;
    return m_v;
}
