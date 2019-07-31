#include "gradientdescent.h"
#include "../WaveFunctions/wavefunction.h"
#include "../sampler.h"
#include "../system.h"
#include <cassert>
#include <iostream>

GradientDescent::GradientDescent(System *system, const double gamma, const double monotonicExp)
    : Optimization(system)
{
    m_numberOfElements = m_system->getNumberOfElements();
    m_maxParameters = m_system->getMaxParameters();
    m_waveFunctionVector = m_system->getWaveFunctionElements();
    m_eta = m_system->getLearningRate();
    m_v = Eigen::MatrixXd::Ones(m_numberOfElements, m_maxParameters);
    m_gamma = gamma;
    m_monotonicExp = monotonicExp;
}

Eigen::MatrixXd GradientDescent::updateParameters()
{
    m_step += 1;
    double monotonic = 1 / pow(m_step, m_monotonicExp);
    m_v = m_gamma * m_v + m_eta * Optimization::getEnergyGradient() * monotonic;
    return m_v;
}
