#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include "sampler.h"
#include "system.h"
#include "Hamiltonians/hamiltonian.h"
#include "WaveFunctions/wavefunction.h"
#include "Optimization/optimization.h"

using std::cout;
using std::endl;


Sampler::Sampler(System* system) {
    m_system                            = system;
    m_numberOfElements                  = m_system->getNumberOfWaveFunctionElements();
    m_maxNumberOfParametersPerElement   = m_system->getMaxNumberOfParametersPerElement();
    m_numberOfMetropolisSteps           = m_system->getNumberOfMetropolisSteps();
    m_equilibriumFraction               = m_system->getEquilibrationFraction();
    m_numberOfStepsAfterEquilibrium     = int((1 - m_system->getEquilibrationFraction()) * m_numberOfMetropolisSteps);
}

void Sampler::sample(const bool acceptedStep, const int stepNumber) {
    if (stepNumber == int(m_numberOfMetropolisSteps * m_equilibriumFraction)) {
        m_acceptenceRatio  = 0;
        m_cumulativeEnergy = 0;
        m_dE    = Eigen::MatrixXd::Zero(m_numberOfElements, m_maxNumberOfParametersPerElement);
        m_dEE   = Eigen::MatrixXd::Zero(m_numberOfElements, m_maxNumberOfParametersPerElement);
        m_SqrdE = 0;
    }

    double EL = m_system->getHamiltonian()->computeLocalEnergy();
    Eigen::MatrixXd gradients = m_system->getOptimization()->getAllImmediateGradients();

    m_cumulativeEnergy  += EL;
    m_dE                += gradients;
    m_dEE               += gradients * EL;
    m_SqrdE             += EL * EL;

    if(acceptedStep) { m_acceptenceRatio += 1; }
}

void Sampler::computeAverages() {
    m_energy   = m_cumulativeEnergy / m_numberOfStepsAfterEquilibrium;
    m_variance = m_SqrdE / m_numberOfStepsAfterEquilibrium - m_energy * m_energy;
}

void Sampler::printOutputToTerminal(const int iter, const int maxIter, const double time) {
    int     ms = m_numberOfMetropolisSteps;
    double  ef = m_equilibriumFraction;

    cout << endl;
    cout << "  -- System info: " << iter+1 << "/" << maxIter << " -- " << endl;
    cout << " Number of particles  : " << m_system->getNumberOfParticles()  << endl;
    cout << " Number of dimensions : " << m_system->getNumberOfDimensions() << endl;
    cout << " Number of Metropolis steps run : 10^" << std::log10(ms) << " (10^" << std::log10(std::round(ms*ef)) << " equilibrium)" << endl;
    cout << endl;
    cout << "  -- Results -- " << endl;
    cout << " Energy           : " << m_energy << endl;
    cout << " Acceptence Ratio : " << double(m_acceptenceRatio)/m_numberOfStepsAfterEquilibrium << endl;
    cout << " Variance         : " << m_variance << endl;
    cout << " STD              : " << sqrt(m_variance) << endl;
    cout << " Time             : " << time << endl;
    cout << endl;
}

std::string generate_filename(const std::string name, const std::string extension) {
    return name + extension;
}

void Sampler::openOutputFiles(const std::string path) {
    std::string energy_filename = generate_filename("energy", ".dat");
    m_energyFile.open(path + energy_filename);                    //Open energy file based on parameters
}

void Sampler::printOutputToFile() {
    m_energyFile << m_energy << endl;;
}

void Sampler::closeOutputFiles() {
    if(m_energyFile.is_open())  { m_energyFile.close(); }
}
