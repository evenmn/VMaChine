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
    m_m = Eigen::MatrixXd::Ones(m_numberOfElements, m_maxParameters);
    m_v = Eigen::MatrixXd::Ones(m_numberOfElements, m_maxParameters);
    m_theta = Eigen::MatrixXd::Ones(m_numberOfElements, m_maxParameters);
}

Eigen::MatrixXd ADAM::updateParameters()
{
    m_step += 1;
    m_g = Optimization::getEnergyGradient();
    m_m = m_beta1 * m_m + (1 - m_beta1) * m_g;
    m_v = m_beta2 * m_v + (1 - m_beta2) * m_g.cwiseAbs2();
    m_mHat = m_m / (1 - pow(m_beta1, m_step));
    m_vHat = m_v / (1 - pow(m_beta2, m_step));
    m_denom = m_vHat.cwiseSqrt() + Eigen::MatrixXd::Ones(m_numberOfElements, m_maxParameters);
    m_theta = m_eta * m_mHat.cwiseProduct(m_denom.cwiseInverse());
    /*for (int i = 0; i < m_numberOfElements; i++) {
        for (int j = 0; j < m_maxParameters; j++) {
            m_theta(i, j) = m_eta * m_mHat(i, j) / (sqrt(m_vHat(i, j)) + m_epsilon);
        }
    }*/
    return m_theta;
}
