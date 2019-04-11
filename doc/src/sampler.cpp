#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <string>
#include "sampler.h"
#include "system.h"
#include "RNG/rng.h"
#include "Hamiltonians/hamiltonian.h"
#include "WaveFunctions/wavefunction.h"
#include "Optimization/optimization.h"
#include "Resampling/AutoBlocking/blocker.h"

using std::cout;
using std::endl;


Sampler::Sampler(System* system) {
    m_system                            = system;
    m_numberOfParticles                 = m_system->getNumberOfParticles();
    m_numberOfDimensions                = m_system->getNumberOfDimensions();
    m_numberOfElements                  = m_system->getNumberOfWaveFunctionElements();
    m_maxNumberOfParametersPerElement   = m_system->getMaxNumberOfParametersPerElement();
    m_numberOfMetropolisSteps           = m_system->getNumberOfMetropolisSteps();
    m_totalNumberOfSteps                = m_system->getTotalNumberOfSteps();
    m_omega                             = m_system->getFrequency();
    m_numberOfBatches                   = m_system->getOptimization()->getNumberOfBatches();
    m_numberOfStepsPerBatch             = int(m_numberOfMetropolisSteps/m_numberOfBatches);
    m_interaction                       = m_system->getInteraction();
}

void Sampler::sample(int numberOfSteps, int equilibriationSteps, const bool acceptedStep, const int stepNumber) {
    if (stepNumber == equilibriationSteps) {
        m_acceptenceRatio           = 0;
        m_cumulativeEnergy          = 0;
        m_cumulativeEnergySqrd      = 0;
        m_cumulativeGradients       = Eigen::MatrixXd::Zero(m_numberOfElements, m_maxNumberOfParametersPerElement);
        m_cumulativeGradientsE      = Eigen::MatrixXd::Zero(m_numberOfElements, m_maxNumberOfParametersPerElement);

        m_equilibriationSteps       = equilibriationSteps;
        m_numberOfSteps             = numberOfSteps;
        m_numberOfStepsPerBatch     = int(m_numberOfSteps/m_numberOfBatches);
    }
    m_instantEnergy    = m_system->getHamiltonian()->computeLocalEnergy();
    m_instantGradients = m_system->getOptimization()->getAllInstantGradients();
    m_cumulativeEnergy      += m_instantEnergy;
    m_cumulativeEnergySqrd  += m_instantEnergy * m_instantEnergy;
    if(stepNumber > (m_numberOfStepsPerBatch * (m_iter%m_numberOfBatches)) && stepNumber < (m_numberOfStepsPerBatch * (m_iter%m_numberOfBatches + 1))) {
        m_cumulativeGradients       += m_instantGradients;
        m_cumulativeGradientsE      += m_instantGradients * m_instantEnergy;
    }
    if(acceptedStep) { m_acceptenceRatio += 1; }
}

void Sampler::computeAverages() {
    m_averageEnergy         = m_cumulativeEnergy / m_numberOfSteps;
    m_averageEnergySqrd     = m_cumulativeEnergySqrd / m_numberOfSteps;
    m_averageGradients      = m_cumulativeGradients / m_numberOfStepsPerBatch;
    m_averageGradientsE     = m_cumulativeGradientsE / m_numberOfStepsPerBatch;
    m_variance              = (m_averageEnergySqrd - m_averageEnergy * m_averageEnergy) / m_numberOfSteps;
}

void Sampler::printOutputToTerminal(const int maxIter, const double time) {
    m_iter += 1;
    cout << endl;
    cout << "  -- System info: " << m_iter << "/" << maxIter << " -- " << endl;
<<<<<<< HEAD
    cout << " Number of particles      : " << m_system->getNumberOfParticles()  << endl;
    cout << " Number of dimensions     : " << m_system->getNumberOfDimensions() << endl;
    cout << " Oscillator frequency     : " << m_omega << endl;
    cout << " # Metropolis steps       : " << m_numberOfSteps + m_equilibriationSteps << " (" << m_numberOfSteps << " equilibration)" << endl;
    cout << " Energy file stored as    : " << m_averageEnergyFileName << endl;
    cout << " Temporary file stored as : " << m_instantEnergyFileName << endl;
=======
    cout << " Number of particles    : " << m_system->getNumberOfParticles()  << endl;
    cout << " Number of dimensions   : " << m_system->getNumberOfDimensions() << endl;
    cout << " Oscillator frequency   : " << m_omega << endl;
    cout << " # Metropolis steps     : " << m_numberOfSteps + m_equilibriationSteps << " (" << m_numberOfSteps << " equilibration)" << endl;
    cout << " Energy file stored as  : " << m_averageEnergyFileName << endl;
    cout << " Instant file stored as : " << m_instantEnergyFileName << endl;
>>>>>>> 66be0ec597973bed31c0d646258f7cfe718f38c3
    cout << endl;
    cout << "  -- Results -- " << endl;
    cout << " Energy            : " << m_averageEnergy << endl;
    cout << " Acceptence Ratio  : " << double(m_acceptenceRatio)/m_numberOfSteps << endl;
    cout << " Variance          : " << m_variance << endl;
    cout << " STD               : " << sqrt(m_variance) << endl;
    cout << " Time              : " << time << endl;
    cout << endl;
}

void Sampler::printFinalOutputToTerminal() {
    std::ifstream infile(m_instantEnergyFileName.c_str());
    std::vector<double> x;
    std::string line;
    while(std::getline(infile,line)) {
        x.push_back(strtod(line.c_str(), nullptr));
    }
    Blocker block(x);

    cout << endl;
    cout << "  ===  Final results:  === " << endl;
    cout << " Energy           : " << m_averageEnergy << endl;
    cout << " Acceptence Ratio : " << double(m_acceptenceRatio)/m_numberOfSteps << endl;
    cout << " Variance         : " << m_variance << endl;
    cout << " STD              : " << sqrt(m_variance) << endl;
    cout << endl;
    cout << " --- Blocking results: ---" << endl;
    printf( " Energy           : %g (with mean sq. err. = %g) \n", block.mean, block.mse_mean);
    printf( " STD              : %g (with mean sq. err. = %g) \n", block.stdErr, block.mse_stdErr);

    if(remove(m_instantEnergyFileName.c_str()) != 0 )
      perror( " Could not remove blocking file" );
    else
      puts( " Removed blocking file" );
}

std::string Sampler::generateFileName(std::string path, std::string name, std::string element, std::string optimization, std::string extension) {
    std::string filename = path;
    filename += "int" + std::to_string(m_interaction) + "/";
    filename += name + "/";
    filename += element + "/";
    filename += std::to_string(m_numberOfDimensions) + "D/";
    filename += std::to_string(m_numberOfParticles) + "P/";
    filename += std::to_string(m_omega) + "w/";
    filename += optimization;
    filename += "_MC" + std::to_string(m_numberOfMetropolisSteps);
    filename += extension;
    return filename;
}

void Sampler::openOutputFiles(const std::string path) {
    // Print average energies to file
    m_averageEnergyFileName = generateFileName(path, "energy", "RBM", "SGD", ".dat");
    m_averageEnergyFile.open(m_averageEnergyFileName);

    // Print instant energies to file
    m_instantEnergyFileName = path + "instant_" + std::to_string(m_system->getRandomNumberGenerator()->nextInt(1e6)) + ".dat";
    m_instantEnergyFile.open(m_instantEnergyFileName);

    // Print onebody densities to file
    if(m_calculateOneBody) {
        std::string oneBodyFileName = generateFileName(path, "onebody", "RBM", "SGD", ".dat");
        m_oneBodyFile.open (oneBodyFileName);
    }
}

void Sampler::printOutputToFile() {
    m_averageEnergyFile << m_averageEnergy << endl;
    if(m_calculateOneBody){
        m_oneBodyFile << m_particlesPerBin << endl;
    }
}

void Sampler::closeOutputFiles() {
    if(m_averageEnergyFile.is_open())      { m_averageEnergyFile.close(); }
    if(m_oneBodyFile.is_open())            { m_oneBodyFile.close(); }
    if(m_instantEnergyFile.is_open())      { m_instantEnergyFile.close(); }
}

void Sampler::printInstantValuesToFile(const Eigen::VectorXd positions) {
    m_instantEnergyFile << m_instantEnergy << endl;  // Write instant energies to file for blocking
    if(m_calculateOneBody) {                         // Calculate onebody densities
        for(int j=0; j<m_numberOfParticles; j++) {
            double dist = 0;
            for(int d=0; d<m_numberOfDimensions; d++) {
                dist += positions(m_numberOfDimensions*j+d) * positions(m_numberOfDimensions*j+d);
            }
            double r = sqrt(dist);      // Distance from particle j to origin
            for(int k=0; k<m_numberOfBins; k++) {
                if(r < m_binLinSpace(k)) {
                    m_particlesPerBin(k) += 1;
                    break;
                }
            }
        }
    }
}
