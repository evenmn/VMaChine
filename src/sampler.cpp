#include <mpi.h>
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
#include "Resampling/blocker.h"

using std::cout;
using std::endl;


Sampler::Sampler(System* system) {
    m_system                            = system;
    m_numberOfProcesses                 = m_system->getNumberOfProcesses();
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
    m_calculateOneBody                  = m_system->getDensity();
    m_printEnergyToFile                 = m_system->getPrintEnergy();
    m_printInstantEnergyToFile          = m_system->getPrintInstantEnergy();
    m_numberOfBins                      = m_system->getNumberOfBins();
    m_maxRadius                         = m_system->getMaxRadius();
    m_radialStep                        = m_maxRadius/m_numberOfBins;
    m_binLinSpace                       = Eigen::VectorXd::LinSpaced(m_numberOfBins, 0, m_maxRadius);
    m_particlesPerBin                   = Eigen::VectorXd::Zero(m_numberOfBins);
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
    m_instantGradients = m_system->getAllInstantGradients();
    m_cumulativeEnergy      += m_instantEnergy;
    m_cumulativeEnergySqrd  += m_instantEnergy * m_instantEnergy;
    if(stepNumber > (m_numberOfStepsPerBatch * (m_iter%m_numberOfBatches)) && stepNumber < (m_numberOfStepsPerBatch * (m_iter%m_numberOfBatches + 1))) {
        m_cumulativeGradients   += m_instantGradients;
        m_cumulativeGradientsE  += m_instantGradients * m_instantEnergy;
    }
    if(acceptedStep) { m_acceptenceRatio += 1; }
}

void Sampler::computeTotals() {
    m_totalCumulativeGradients       = Eigen::MatrixXd::Zero(m_numberOfElements, m_maxNumberOfParametersPerElement);
    m_totalCumulativeGradientsE      = Eigen::MatrixXd::Zero(m_numberOfElements, m_maxNumberOfParametersPerElement);
    MPI_Reduce(&m_cumulativeEnergy,     &m_totalCumulativeEnergy,     1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&m_cumulativeEnergySqrd, &m_totalCumulativeEnergySqrd, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    for(int i=0; i<m_numberOfElements; i++) {
        for(int j=0; j<m_maxNumberOfParametersPerElement; j++) {
            MPI_Reduce(&m_cumulativeGradients(i,j),  &m_totalCumulativeGradients(i,j),  1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(&m_cumulativeGradientsE(i,j), &m_totalCumulativeGradientsE(i,j), 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }
    }
}

void Sampler::computeAverages() {
    int    totalNumberOfMCSamples     = m_numberOfSteps * m_numberOfProcesses;
    m_averageEnergy         = m_totalCumulativeEnergy     / totalNumberOfMCSamples;
    m_averageEnergySqrd     = m_totalCumulativeEnergySqrd / totalNumberOfMCSamples;
    m_averageGradients      = m_totalCumulativeGradients  / m_numberOfStepsPerBatch;
    m_averageGradientsE     = m_totalCumulativeGradientsE / m_numberOfStepsPerBatch;
    m_variance              = (m_averageEnergySqrd - m_averageEnergy * m_averageEnergy) / totalNumberOfMCSamples;
}

void Sampler::printOutputToTerminal(const int maxIter, const double time) {
    m_iter += 1;
    cout << endl;
    cout << "  -- System info: " << m_iter << "/" << maxIter << " -- " << endl;
    cout << " Number of particles    : " << m_system->getNumberOfParticles()  << endl;
    cout << " Number of dimensions   : " << m_system->getNumberOfDimensions() << endl;
    cout << " Number of processes    : " << m_system->getNumberOfProcesses()  << endl;
    cout << " Number of parameters   : " << m_system->getTotalNumberOfParameters() << endl;
    cout << " Oscillator frequency   : " << m_omega << endl;
    cout << " # Metropolis steps     : " << (m_numberOfSteps + m_equilibriationSteps) * m_numberOfProcesses << " ("
                                         << m_numberOfSteps * m_numberOfProcesses << " equilibration)" << endl;
    cout << " Energy file stored as  : " << m_averageEnergyFileName << endl;
    cout << " Instant file stored as : " << m_instantEnergyFileName << endl;
    cout << endl;
    cout << "  -- Results -- " << endl;
    cout << " Energy            : " << m_averageEnergy << endl;
    cout << " Acceptence Ratio  : " << double(m_acceptenceRatio)/m_numberOfSteps << endl;
    cout << " Variance          : " << m_variance << endl;
    cout << " STD               : " << sqrt(m_variance) << endl;
    cout << " CPU Time          : " << time << endl;
    cout << endl;
}

void Sampler::printFinalOutputToTerminal(int instantNumber, std::string path) {
    cout << endl;
    cout << "  ===  Final results:  === " << endl;
    cout << " Energy           : " << m_averageEnergy << endl;
    cout << " Acceptence Ratio : " << double(m_acceptenceRatio)/m_numberOfSteps << endl;
    cout << " Variance         : " << m_variance << endl;
    cout << " STD              : " << sqrt(m_variance) << endl;
    cout << endl;

    if(m_printInstantEnergyToFile) {
        // Append all instant files into one file
        std::ofstream ofile(m_instantEnergyFileName.c_str(), std::ios::out | std::ios::app);
        for(int i=1; i<m_numberOfProcesses; i++) {
            std::string name = path + "instant_" + std::to_string(instantNumber) + "_" + std::to_string(i) + ".dat";
            std::ifstream ifile(name.c_str(), std::ios::in);
            if (!ifile.is_open()) {
                perror ( "File not found" );
            }
            else {
                ofile << ifile.rdbuf();
            }
            if(remove(name.c_str()) != 0)
                perror( " Could not remove blocking file" );
        }
        std::ifstream infile(m_instantEnergyFileName.c_str());
        std::vector<double> x;
        std::string line;
        while(std::getline(infile,line)) {
            x.push_back(strtod(line.c_str(), nullptr));
        }
        Blocker block(x);

        cout << " --- Blocking results: ---" << endl;
        printf( " Energy           : %g (with mean sq. err. = %g) \n", block.mean, block.mse_mean);
        printf( " STD              : %g (with mean sq. err. = %g) \n", block.stdErr, block.mse_stdErr);

        if(remove(m_instantEnergyFileName.c_str()) != 0)
            perror( " Could not remove blocking file" );
        else
            puts( " Successfully removed blocking file" );

        // Merge all one-body files into one file
        std::ifstream infile1;
        std::string mainName = generateFileName(path, "onebody", "SGD", "_" + std::to_string(0) + ".dat");
        infile1.open(mainName.c_str());
        std::ofstream outfile;
        outfile.open("file.dat");
        for(int i=1; i<m_numberOfProcesses; i++) {
            std::string name = generateFileName(path, "onebody", "SGD", "_" + std::to_string(i) + ".dat");
            std::ifstream infile2;
            infile2.open(name.c_str());
            if (!infile1.is_open() || !infile2.is_open()) {
                cout << "file not found";
            }
            else {
                int value, value2;
                while (infile1 >> value && infile2 >> value2) {
                    outfile << double(value) + double(value2) << endl;
                }
            }
            std::rename("file.dat", mainName.c_str());

        }
    }
}

std::string Sampler::generateFileName(std::string path, std::string name, std::string optimization, std::string extension) {
    std::string filename = path;
    filename += "int" + std::to_string(m_interaction) + "/";
    filename += name + "/";
    filename += m_system->getAllLabels() + "/";
    filename += std::to_string(m_numberOfDimensions) + "D/";
    filename += std::to_string(m_numberOfParticles) + "P/";
    filename += std::to_string(m_omega) + "w/";
    filename += optimization;
    filename += "_MC" + std::to_string(m_numberOfMetropolisSteps * m_numberOfProcesses);
    filename += extension;
    return filename;
}

void Sampler::openOutputFiles(const std::string path, int instantNumber, int myRank) {
    // Print average energies to file
    if(m_printEnergyToFile) {
        m_averageEnergyFileName = generateFileName(path, "energy", "SGD", ".dat");
        m_averageEnergyFile.open(m_averageEnergyFileName);
    }

    if(m_printInstantEnergyToFile) {
        m_instantEnergyFileName = path + "instant_" + std::to_string(instantNumber) + "_" + std::to_string(myRank) + ".dat";
        m_instantEnergyFile.open(m_instantEnergyFileName);
    }
    if(m_calculateOneBody) {
        std::string oneBodyFileName = generateFileName(path, "onebody", "SGD", "_" + std::to_string(myRank) + ".dat");
        m_oneBodyFile.open (oneBodyFileName);
    }
}

void Sampler::printOutputToFile(int myRank) {
    if(m_printEnergyToFile && myRank == 0) {
        m_averageEnergyFile << m_averageEnergy << endl;
    }
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
    if(m_printInstantEnergyToFile) {
        m_instantEnergyFile << m_instantEnergy << endl;  // Write instant energies to file for blocking
    }
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
