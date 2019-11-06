#include "optimization.h"
#include "../sampler.h"
#include "../system.h"

Optimization::Optimization(System *system)
{
    m_system = system;
}

arma::mat Optimization::getEnergyGradient()
{
    double averageEnergy = m_system->getSampler()->getAverageEnergy();
    arma::mat averageGradients = m_system->getSampler()->getAverageGradients();
    arma::mat averageGradientsE = m_system->getSampler()->getAverageGradientsE();
    return 2 * (averageGradientsE - averageEnergy * averageGradients);
}

Optimization::~Optimization() {}
