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
#include "block/c++/blocker.h"

using std::cout;
using std::endl;


Sampler::Sampler(System* system) {
    m_system                            = system;
    m_numberOfProcesses                 = m_system->getNumberOfProcesses();
    m_numberOfParticles                 = m_system->getNumberOfParticles();
    m_numberOfDimensions                = m_system->getNumberOfDimensions();
    m_numberOfElements                  = m_system->getNumberOfWaveFunctionElements();
    m_maxNumberOfParametersPerElement   = m_system->getMaxNumberOfParametersPerElement();

    m_totalNumberOfStepsWOEqui          = m_system->getTotalNumberOfStepsWOEqui();
    m_totalNumberOfStepsWEqui           = m_system->getTotalNumberOfStepsWEqui();
    m_numberOfStepsWOEqui               = m_system->getNumberOfStepsWOEqui();
    m_initialTotalNumberOfStepsWOEqui   = m_system->getInitialTotalNumberOfStepsWOEqui();
    m_numberOfEquilibriationSteps       = m_system->getnumberOfEquilibriationSteps();

    m_omega                             = m_system->getFrequency();
    m_numberOfBatches                   = m_system->getOptimization()->getNumberOfBatches();
    m_numberOfStepsPerBatch             = unsigned(m_numberOfStepsWOEqui/m_numberOfBatches);
    m_interaction                       = m_system->getInteraction();
    m_computeOneBodyDensity             = m_system->getDensity();
    m_computeTwoBodyDensity             = m_system->computeTwoBodyDensity();
    m_printEnergyToFile                 = m_system->getPrintEnergy();
    m_printInstantEnergyToFile          = m_system->getPrintInstantEnergy();
    m_numberOfBins                      = m_system->getNumberOfBins();
    m_maxRadius                         = m_system->getMaxRadius();
    m_rank                              = m_system->getRank();
    m_path                              = m_system->getPath();
    m_radialStep                        = m_maxRadius/m_numberOfBins;
    m_binLinSpace                       = Eigen::VectorXd::LinSpaced(m_numberOfBins, 0, m_maxRadius);
    m_particlesPerBin                   = Eigen::VectorXi::Zero(m_numberOfBins);
    m_particlesPerBinPairwise           = Eigen::MatrixXi::Zero(m_numberOfBins, m_numberOfBins);
}

void Sampler::sample(const bool acceptedStep, const unsigned long stepNumber) {
    if (stepNumber == m_numberOfEquilibriationSteps) {
        m_acceptence                = 0;
        m_cumulativeEnergy          = 0;
        m_cumulativeEnergySqrd      = 0;
        m_cumulativeGradients       = Eigen::MatrixXd::Zero(m_numberOfElements, m_maxNumberOfParametersPerElement);
        m_cumulativeGradientsE      = Eigen::MatrixXd::Zero(m_numberOfElements, m_maxNumberOfParametersPerElement);
    }
    m_instantEnergy    = m_system->getHamiltonian()->computeLocalEnergy();
    m_instantGradients = m_system->getAllInstantGradients();
    m_cumulativeEnergy      += m_instantEnergy;
    m_cumulativeEnergySqrd  += m_instantEnergy * m_instantEnergy;
    if(stepNumber > (m_numberOfStepsPerBatch * (m_iter%m_numberOfBatches)) && stepNumber < (m_numberOfStepsPerBatch * (m_iter%m_numberOfBatches + 1))) {
        m_cumulativeGradients   += m_instantGradients;
        m_cumulativeGradientsE  += m_instantGradients * m_instantEnergy;
    }
    if(acceptedStep) { m_acceptence += 1; }
}

void Sampler::computeTotals() {
    m_totalCumulativeGradients       = Eigen::MatrixXd::Zero(m_numberOfElements, m_maxNumberOfParametersPerElement);
    m_totalCumulativeGradientsE      = Eigen::MatrixXd::Zero(m_numberOfElements, m_maxNumberOfParametersPerElement);
    MPI_Reduce(&m_acceptence,           &m_totalAcceptence,           1, MPI_INT,    MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&m_cumulativeEnergy,     &m_totalCumulativeEnergy,     1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&m_cumulativeEnergySqrd, &m_totalCumulativeEnergySqrd, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    for(unsigned int i=0; i<m_numberOfElements; i++) {
        for(unsigned int j=0; j<m_maxNumberOfParametersPerElement; j++) {
            MPI_Reduce(&m_cumulativeGradients(i,j),  &m_totalCumulativeGradients(i,j),  1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(&m_cumulativeGradientsE(i,j), &m_totalCumulativeGradientsE(i,j), 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }
    }
}

void Sampler::setNumberOfSteps(const unsigned long numberOfStepsWOEqui, const unsigned long totalNumberOfStepsWOEqui, const unsigned long totalNumberOfStepsWEqui) {
    m_numberOfStepsWOEqui      = numberOfStepsWOEqui;
    m_totalNumberOfStepsWOEqui = totalNumberOfStepsWOEqui;
    m_totalNumberOfStepsWEqui  = totalNumberOfStepsWEqui;
    m_numberOfStepsPerBatch    = unsigned(m_totalNumberOfStepsWOEqui/m_numberOfBatches);
}

void Sampler::computeAverages() {
    m_averageEnergy         = m_totalCumulativeEnergy     / m_totalNumberOfStepsWOEqui;
    m_averageEnergySqrd     = m_totalCumulativeEnergySqrd / m_totalNumberOfStepsWOEqui;
    m_averageGradients      = m_totalCumulativeGradients  / m_numberOfStepsPerBatch;
    m_averageGradientsE     = m_totalCumulativeGradientsE / m_numberOfStepsPerBatch;
    m_variance              = (m_averageEnergySqrd - m_averageEnergy * m_averageEnergy) / m_totalNumberOfStepsWOEqui;
}

void Sampler::printOutputToTerminal(const unsigned int maxIter, const double time) {
    m_iter += 1;
    cout << endl;
    cout << "  -- System info: " << " -- " << endl;
    cout << " Iteration progress      : " << m_iter << "/" << maxIter << endl;
    cout << " Number of particles     : " << m_system->getNumberOfParticles()  << endl;
    cout << " Number of dimensions    : " << m_system->getNumberOfDimensions() << endl;
    cout << " Number of processes     : " << m_system->getNumberOfProcesses()  << endl;
    cout << " Number of parameters    : " << m_system->getTotalNumberOfParameters() << endl;
    cout << " Oscillator frequency    : " << m_omega << endl;
    cout << " # Metropolis steps      : " << m_totalNumberOfStepsWEqui  << " ("
                                          << m_totalNumberOfStepsWOEqui << " equilibration)" << endl;
    cout << " Energy file stored as   : " << m_averageEnergyFileName << endl;
    cout << " Blocking file stored as : " << m_instantEnergyFileName << endl;
    cout << endl;
    cout << "  -- Results -- " << endl;
    cout << " Energy            : " << m_averageEnergy << endl;
    cout << " Acceptence Ratio  : " << double(m_totalAcceptence)/m_totalNumberOfStepsWOEqui << endl;
    cout << " Variance          : " << m_variance << endl;
    cout << " STD               : " << sqrt(m_variance) << endl;
    cout << " CPU Time          : " << time << endl;
    cout << endl;
}

void Sampler::printFinalOutputToTerminal() {
    cout << endl;
    cout << "  ===  Final results:  === " << endl;
    cout << " Energy           : " << m_averageEnergy << endl;
    cout << " Acceptence Ratio : " << double(m_totalAcceptence)/m_totalNumberOfStepsWOEqui << endl;
    cout << " Variance         : " << m_variance << endl;
    cout << " STD              : " << sqrt(m_variance) << endl;
    cout << endl;

    if(m_printInstantEnergyToFile) {
        appendInstantFiles();

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
    }
}

void Sampler::appendInstantFiles() {
    std::ofstream outfile(m_instantEnergyFileName.c_str(), std::ios::out | std::ios::app);
    for(int i=1; i<m_numberOfProcesses; i++) {
        std::string name = m_path + "instant_" + std::to_string(m_instantNumber) + "_" + std::to_string(i) + ".dat";
        std::ifstream infile(name.c_str(), std::ios::in);
        if (!infile.is_open()) {
            perror ( "File not found" );
        }
        else {
            outfile << infile.rdbuf();
        }
        if(remove(name.c_str()) != 0)
            perror( " Could not remove blocking file" );
    }
}

std::string Sampler::generateFileName(std::string path, std::string name, std::string optimization, std::string extension) {
    std::string fileName = path;
    fileName += "int" + std::to_string(m_interaction) + "/";
    fileName += name + "/";
    fileName += m_system->getAllLabels() + "/";
    fileName += std::to_string(m_numberOfDimensions) + "D/";
    fileName += std::to_string(m_numberOfParticles) + "P/";
    fileName += std::to_string(m_omega) + "w/";
    fileName += optimization;
    fileName += "_MC" + std::to_string(m_initialTotalNumberOfStepsWOEqui);
    fileName += extension;
    return fileName;
}

void Sampler::openOutputFiles() {
    if(m_rank == 0) {
        m_instantNumber = m_system->getRandomNumberGenerator()->nextInt(unsigned(1e6));
    }
    MPI_Bcast(&m_instantNumber, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Print average energies to file
    if(m_printEnergyToFile && m_rank == 0) {
        m_averageEnergyFileName = generateFileName(m_path, "energy", "SGD", ".dat");
        m_averageEnergyFile.open(m_averageEnergyFileName);
    }

    if(m_printInstantEnergyToFile) {
        m_instantEnergyFileName = m_path + "instant_" + std::to_string(m_instantNumber) + "_" + std::to_string(m_rank) + ".dat";
        m_instantEnergyFile.open(m_instantEnergyFileName);
    }
    if(m_computeOneBodyDensity && m_rank == 0) {
        std::string oneBodyFileName = generateFileName(m_path, "onebody", "SGD", ".dat");
        m_oneBodyFile.open (oneBodyFileName);
    }
    if(m_computeTwoBodyDensity && m_rank == 0) {
        std::string twoBodyFileName = generateFileName(m_path, "twobody", "SGD", ".dat");
        m_twoBodyFile.open (twoBodyFileName);
    }
}

void Sampler::printEnergyToFile() {
    if(m_printEnergyToFile && m_rank == 0) {
        m_averageEnergyFile << m_averageEnergy << endl;
    }
}

void Sampler::printOneBodyDensityToFile() {
    if(m_computeOneBodyDensity){
        m_totalParticlesPerBin = Eigen::VectorXi::Zero(m_numberOfBins);
        //MPI_Gather(m_particlesPerBin.data(), m_numberOfBins, MPI_DOUBLE, m_totalParticlesPerBin.data(), m_numberOfBins, MPI_DOUBLE, MPI_COMM_WORLD);
        for(unsigned int i=0; i<m_numberOfBins; i++) {
            MPI_Reduce(&m_particlesPerBin(i), &m_totalParticlesPerBin(i), 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        }
        if(m_rank == 0) {
            m_oneBodyFile << m_totalParticlesPerBin << endl;
        }
    }
}

void Sampler::printTwoBodyDensityToFile() {
    if(m_computeTwoBodyDensity){
        m_totalParticlesPerBinPairwise = Eigen::MatrixXi::Zero(m_numberOfBins, m_numberOfBins);
        for(unsigned int i=0; i<m_numberOfBins; i++) {
            for(unsigned int j=0; j<m_numberOfBins; j++) {
                MPI_Reduce(&m_particlesPerBinPairwise(i,j), &m_totalParticlesPerBinPairwise(i,j), 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            }
        }
        if(m_rank == 0) {
            m_twoBodyFile << m_totalParticlesPerBinPairwise << endl;
        }
    }
}

void Sampler::closeOutputFiles() {
    if(m_averageEnergyFile.is_open())      { m_averageEnergyFile.close(); }
    if(m_oneBodyFile.is_open())            { m_oneBodyFile.close(); }
    if(m_twoBodyFile.is_open())            { m_twoBodyFile.close(); }
    if(m_instantEnergyFile.is_open())      { m_instantEnergyFile.close(); }
}

void Sampler::printInstantValuesToFile() {
    if(m_printInstantEnergyToFile) {
        m_instantEnergyFile << m_instantEnergy << endl;  // Write instant energies to file for blocking
    }
}

void Sampler::computeOneBodyDensity(const Eigen::VectorXd positions) {
    if(m_computeOneBodyDensity) {                         // Calculate onebody densities
        for(unsigned int particle=0; particle<m_numberOfParticles; particle++) {
            double dist = 0;
            unsigned int positionIndex = m_numberOfDimensions * particle;
            for(unsigned short d=0; d<m_numberOfDimensions; d++) {
                double position = positions(positionIndex+d);
                dist += position * position;
            }
            double r = sqrt(dist);      // Distance from particle to origin
            for(unsigned int k=0; k<m_numberOfBins; k++) {
                if(r < m_binLinSpace(k)) {
                    m_particlesPerBin(k) += 1;
                    break;
                }
            }
        }
    }
}

void Sampler::computeTwoBodyDensity(const Eigen::VectorXd positions) {
    if(m_computeTwoBodyDensity) {                         // Calculate twobody densities
        for(unsigned int particle1=0; particle1<m_numberOfParticles; particle1++) {
            double dist1 = 0;
            unsigned int position1Index = m_numberOfDimensions * particle1;
            for(unsigned short d=0; d<m_numberOfDimensions; d++) {
                double position1 = positions(position1Index+d);
                dist1 += position1 * position1;
            }
            double r1 = sqrt(dist1);      // Distance from particle 1 to origin
            int counter1 = 0;
            while(m_binLinSpace(counter1) < r1) {
                counter1 += 1;
            }
            for(unsigned int particle2=particle1+1; particle2<m_numberOfParticles; particle2++) {
                double dist2 = 0;
                unsigned int position2Index = m_numberOfDimensions * particle2;
                for(unsigned short d=0; d<m_numberOfDimensions; d++) {
                    double position2 = positions(position2Index+d);
                    dist2 += position2 * position2;
                }
                double r2 = sqrt(dist2);      // Distance from particle 2 to origin
                int counter2 = 0;
                while(m_binLinSpace(counter2) < r2) {
                    counter2 += 1;
                }
                m_particlesPerBinPairwise(counter1, counter2) += 1;
            }
        }
    }
}
