#include "adam.h"
#include "../WaveFunctions/wavefunction.h"
#include "../sampler.h"
#include "../system.h"
#include <cassert>
#include <iostream>

ADAM::ADAM(System *system)
    : Optimization(system)
{}

void ADAM::initialize()
{
    m_numberOfElements = m_system->getNumberOfElements();
    m_maxParameters = m_system->getMaxParameters();
    m_eta = m_system->getLearningRate();
    //arma::uvec = (m_numberOfElements, m_maxParameters);
    m_m.ones(m_numberOfElements, m_maxParameters);
    m_v.ones(m_numberOfElements, m_maxParameters);
    m_theta.ones(m_numberOfElements, m_maxParameters);
}

arma::mat ADAM::updateParameters()
{
    m_step += 1;
    m_g = Optimization::getEnergyGradient();
    m_m = m_beta1 * m_m + (1 - m_beta1) * m_g;
    m_v = m_beta2 * m_v + (1 - m_beta2) * arma::square(m_g);
    m_mHat = m_m / (1 - pow(m_beta1, m_step));
    m_vHat = m_v / (1 - pow(m_beta2, m_step));
    m_denom = arma::sqrt(m_vHat) + arma::ones(size(m_vHat));
    m_theta = m_eta * m_mHat % arma::pow(m_denom, -1);
    /*for (arma::uword i = 0; i < m_numberOfElements; i++) {
        for (arma::uword j = 0; j < m_maxParameters; j++) {
            m_theta(i, j) = m_eta * m_mHat(i, j) / (sqrt(m_vHat(i, j)) + m_epsilon);
        }
    }*/
    return m_theta;
}
