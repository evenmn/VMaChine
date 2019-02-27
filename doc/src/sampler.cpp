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
    m_numberOfParticles                 = m_system->getNumberOfParticles();
    m_numberOfDimensions                = m_system->getNumberOfDimensions();
    m_numberOfElements                  = m_system->getNumberOfWaveFunctionElements();
    m_maxNumberOfParametersPerElement   = m_system->getMaxNumberOfParametersPerElement();
    m_numberOfMetropolisSteps           = m_system->getNumberOfMetropolisSteps();
    m_equilibriumFraction               = m_system->getEquilibrationFraction();
    m_totalNumberOfSteps                = m_system->getTotalNumberOfSteps();
    m_numberOfMetropolisSteps           = m_system->getNumberOfMetropolisSteps();
}

void Sampler::sample(const bool acceptedStep, const int stepNumber) {
    if (stepNumber == (m_totalNumberOfSteps - m_numberOfMetropolisSteps)) {
        m_acceptenceRatio  = 0;
        m_cumulativeEnergy = 0;
        m_dE    = Eigen::MatrixXd::Zero(m_numberOfElements, m_maxNumberOfParametersPerElement);
        m_dEE   = Eigen::MatrixXd::Zero(m_numberOfElements, m_maxNumberOfParametersPerElement);
        m_SqrdE = 0;
    }

    m_instantEnergy = m_system->getHamiltonian()->computeLocalEnergy();
    Eigen::MatrixXd gradients = m_system->getOptimization()->getAllImmediateGradients();

    m_cumulativeEnergy  += m_instantEnergy;
    m_dE                += gradients;
    m_dEE               += gradients * m_instantEnergy;
    m_SqrdE             += m_instantEnergy * m_instantEnergy;

    if(acceptedStep) { m_acceptenceRatio += 1; }
}

void Sampler::computeAverages() {
    m_energy   = m_cumulativeEnergy / m_numberOfMetropolisSteps;
    m_variance = m_SqrdE / m_numberOfMetropolisSteps - m_energy * m_energy;
}

void Sampler::printOutputToTerminal(const int iter, const int maxIter, const double time) {
    cout << endl;
    cout << "  -- System info: " << iter+1 << "/" << maxIter << " -- " << endl;
    cout << " Number of particles  : " << m_system->getNumberOfParticles()  << endl;
    cout << " Number of dimensions : " << m_system->getNumberOfDimensions() << endl;
    cout << " Number of Metropolis steps run : " << m_totalNumberOfSteps << " (" << m_numberOfMetropolisSteps << " equilibrium)" << endl;
    cout << endl;
    cout << "  -- Results -- " << endl;
    cout << " Energy           : " << m_energy << endl;
    cout << " Acceptence Ratio : " << double(m_acceptenceRatio)/m_numberOfMetropolisSteps << endl;
    cout << " Variance         : " << m_variance << endl;
    cout << " STD              : " << sqrt(m_variance) << endl;
    cout << " Time             : " << time << endl;
    cout << endl;
}

std::string generateFileName(const std::string name, const std::string extension) {
    return name + extension;
}

void Sampler::openOutputFiles(const std::string path) {
    // Print average energies to file
    std::string energyFileName = generateFileName("energy", ".dat");
    m_energyFile.open(path + energyFileName);

    // Print cumulative energies to file
    std::string cumulativeFileName = generateFileName("cumulative", ".dat");
    m_cumulativeFile.open(path + cumulativeFileName);

    // Print onebody densities to file
    if(m_calculateOneBody) {
        std::string oneBodyFileName = generateFileName("OB", ".dat");
        m_oneBodyFile.open (path + oneBodyFileName);
    }
}

void Sampler::printOutputToFile() {
    m_energyFile << m_energy << endl;
    if(m_calculateOneBody){
        m_oneBodyFile << m_particlesPerBin << endl;
    }
}

void Sampler::closeOutputFiles() {
    if(m_energyFile.is_open())      { m_energyFile.close(); }
    if(m_oneBodyFile.is_open())     { m_oneBodyFile.close(); }
    if(m_cumulativeFile.is_open())  { m_cumulativeFile.close(); }
}

void Sampler::printImmediatelyToFile(const Eigen::VectorXd positions) {
    // Write cumulative energies to file for blocking
    m_cumulativeFile << m_instantEnergy << endl;

    // Calculate onebody densities
    if(m_calculateOneBody) {
        for(int j=0; j<m_numberOfParticles; j++) {
            double dist = 0;
            for(int d=0; d<m_numberOfDimensions; d++) {
                dist += positions(m_numberOfDimensions*j+d) * positions(m_numberOfDimensions*j+d);
            }
            double r = sqrt(dist);      //Distance from particle j to origin
            for(int k=0; k<m_numberOfBins; k++) {
                if(r < m_binLinSpace(k)) {
                    m_particlesPerBin(k) += 1;
                    break;
                }
            }
        }
    }
}
