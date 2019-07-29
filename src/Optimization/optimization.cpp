#include "optimization.h"
#include "../sampler.h"
#include "../system.h"

Optimization::Optimization(System *system)
{
    m_system = system;
}

Eigen::MatrixXd Optimization::getEnergyGradient()
{
    double averageEnergy = m_system->getSampler()->getAverageEnergy();
    Eigen::MatrixXd averageGradients = m_system->getSampler()->getAverageGradients();
    Eigen::MatrixXd averageGradientsE = m_system->getSampler()->getAverageGradientsE();
    return 2 * (averageGradientsE - averageEnergy * averageGradients);
}

Optimization::~Optimization() {}
